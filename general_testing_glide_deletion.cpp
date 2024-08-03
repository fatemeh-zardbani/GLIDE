#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>

#include "GLIDE_GSM.h"

#ifdef fo8
    const int fo = 8;
#elif fo16
    const int fo = 16;
#elif fo32
    const int fo = 32;
#endif


#ifdef maxt64
    const int maxt = 64;
#elif maxt256
    const int maxt = 256;
#elif maxt512
    const int maxt = 512;
#elif maxt1024
    const int maxt = 1024;
#elif maxt2048
    const int maxt = 2048;
#endif


#ifdef mint32
    const int mint = 32;
#elif mint48
    const int mint = 48;
#elif mint128
    const int mint = 128;
#elif mint192
    const int mint = 192;
#elif mint256
    const int mint = 256;
#elif mint384
    const int mint = 384;
#elif mint512
    const int mint = 512;
#elif mint768
    const int mint = 768;
#elif mint1024
    const int mint = 1024;
#elif mint1536
    const int mint = 1536;
#endif

using namespace std;


typedef int ValueType;
typedef float ELEMTYPE;
typedef double ELEMTYPEREAL;

struct Rect
{
    Rect()  {}

    Rect(ELEMTYPE a_min[NUMDIMS], ELEMTYPE a_max[NUMDIMS]){
        for(int i = 0; i < NUMDIMS; i++) m_min[i] = a_min[i];
        for(int i = 0; i < NUMDIMS; i++) m_max[i] = a_max[i];
    }
    ELEMTYPE m_min[NUMDIMS];
    ELEMTYPE m_max[NUMDIMS];
};


bool Overlap(Rect *a_rectA, Rect *a_rectB);

bool Overlap(Rect *a_rectA, Rect *a_rectB) {
            ASSERT(a_rectA && a_rectB);

    for (int index = 0; index < 2; ++index) {
        if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
            a_rectB->m_min[index] > a_rectA->m_max[index]) {
            return false;
        }
    }
    return true;
}



double timing(){
    static struct timeval t1, t2;
    gettimeofday(&t2,NULL);
    double ret = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) * 1e-6;
    t1 = t2;
    return ret;
}


int main(int argc, char **argv){
    // auto seed = time(NULL);
    auto seed = 0;
    srand(seed);

    string data_file_name = argv[1];
    int query_size = stoi(argv[2]);
    string query_file_name = argv[3];
    string time_file_name = argv[4];
    int ratio = stoi(argv[5]);
    int insert_count = stoi(argv[6]);
    int delete_count = stoi(argv[7]);
    string delete_file_name = argv[8];
    
    std::ifstream query_file(query_file_name.c_str());

    #ifdef stats
        ofstream stats_file;
        string stats_file_name = argv[8];
        stats_file.open(stats_file_name.c_str());   
    #endif

    ofstream times_file;
    times_file.open(time_file_name.c_str());

    
    typedef RTree<ValueType, fo, fo/2, maxt, mint> MyTree;
    MyTree tree(data_file_name);

    ELEMTYPE min[NUMDIMS]; ELEMTYPE max[NUMDIMS];
    int this_id;

    query_file.clear();
    query_file.seekg(0, ios::beg);
    Rect queries[query_size];

    vector<Rect> data(DATA_COUNT);
    vector<int> data_ids(DATA_COUNT);
    std::ifstream data_file(data_file_name.c_str());

    for(int i = 0; i < DATA_COUNT; i++){
        data[i] = Rect();
        data_ids[i] = i;
        for(int j = 0; j < NUMDIMS; j++) data_file >> data[i].m_min[j];
        for(int j = 0; j < NUMDIMS; j++) data_file >> data[i].m_max[j];
    }
    

    for(int i = 0; i < query_size; i++){
        queries[i] = Rect();
        for(int j = 0; j < NUMDIMS; j++) query_file >> queries[i].m_min[j];
        for(int j = 0; j < NUMDIMS; j++) query_file >> queries[i].m_max[j];
       
    }

    int how_many_deleted_till_now = 0;
    std::ifstream delete_file(delete_file_name.c_str());
    vector<Rect> tbd_rects;
    vector<int> tbd_ids;
    int trash_id;
    for(int i = 0; i < delete_count*(query_size/ratio - 1); i++){

        for(int j = 0; j < NUMDIMS; j++) delete_file >> min[j];
        for(int j = 0; j < NUMDIMS; j++) delete_file >> max[j];
        delete_file >> trash_id;
        tbd_rects.push_back(Rect(min, max));
        tbd_ids.push_back(trash_id);
    }
    int found_count; int sum_ids; int total_count;
    int lin_count;
    clock_t q_time = 0;
    int rightestright = 0;
    int last_id = DATA_COUNT -1;
    int item_index_to_be_deleted_choice = 0;
    clock_t insert_time, start_insert_time, delete_time, start_delete_time;
    clock_t start_query_time;
    
    for(int i = 0; i < query_size; i++){
        cerr << "QUERY " << i << "...\n";

        insert_time = 0.0;
        delete_time = 0.0;
        if(i % ratio == 0 && i != 0){
            for(int o =0; o < insert_count; o++){
                for(int j = 0; j < NUMDIMS; j++) data_file >> min[j];
                for(int j = 0; j < NUMDIMS; j++) data_file >> max[j];
                this_id = ++last_id;
                start_insert_time = clock();
                tree.Insert_gradual(min, max, this_id);
                insert_time += clock() - start_insert_time;
                data.push_back(Rect(min, max));
                data_ids.push_back(this_id);
            }
            
            #ifdef delete_timebreak
                spare_scan_time = 0.0;
                leaf_scan_time = 0.0;
            #endif
            for(int o=0; o < delete_count; o++){
                start_delete_time = clock();
                tree.Delete(tbd_rects[how_many_deleted_till_now].m_min, tbd_rects[how_many_deleted_till_now].m_max, tbd_ids[how_many_deleted_till_now]);
                delete_time += clock() - start_delete_time;
                how_many_deleted_till_now++;
            }
            #ifdef delete_timebreak
                cerr << (double)spare_scan_time/CLOCKS_PER_SEC << " " << (double)leaf_scan_time/CLOCKS_PER_SEC << endl;
            #endif
        }
        start_query_time = clock();
        found_count = tree.QueryAdaptive(queries[i].m_min, queries[i].m_max);
        #ifdef stats
            int rightestright = tree.getRightestRight();
            int height = tree.TreeHeight();
            stats_file << "i " << rightestright << " " << count_total_holes << " " << height << " " << count_internal_nodes << " " << count_regular_leaves << " " << count_irregular_leaves << endl;
        #endif
        q_time = clock() - start_query_time;
        
        times_file << found_count << " " << (double)insert_time/CLOCKS_PER_SEC << " " << (double)delete_time/CLOCKS_PER_SEC << " " << (double)q_time/CLOCKS_PER_SEC <<  "\n";

	    cout << found_count << endl;
    }
    query_file.close();
    cout << "DONE!!\n";
    return 0;
}

