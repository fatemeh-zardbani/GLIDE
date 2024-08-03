#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>

#ifdef CQRIPPLE
    #include "GLIDE_CQripple.h"
#elif CQ
    #include "GLIDE_CQ.h"
#elif CQM
    #include "GLIDE_CQM.h"
#elif CS
    #include "GLIDE_CS.h"
#elif CSM
    #include "GLIDE_CSM.h"
#elif GQ
    #include "GLIDE_GQ.h"
#elif GQM
    #include "GLIDE_GQM.h"
#elif GS
    #include "GLIDE_GS.h"
#elif GSM
    #include "GLIDE_GSM.h"
#elif GSQ
    #include "GLIDE_GSQ.h"
#elif GSS
    #include "GLIDE_GSS.h"
#endif



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
    
    std::ifstream query_file(query_file_name.c_str());

    #ifdef stats
        ofstream stats_file;
        string stats_file_name = argv[7];
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
    std::ifstream data_file(data_file_name.c_str());

    for(int i = 0; i < DATA_COUNT; i++){
        data[i] = Rect();
        for(int j = 0; j < NUMDIMS; j++) data_file >> data[i].m_min[j];
        for(int j = 0; j < NUMDIMS; j++) data_file >> data[i].m_max[j];
    }
    

    for(int i = 0; i < query_size; i++){
        queries[i] = Rect();
        for(int j = 0; j < NUMDIMS; j++) query_file >> queries[i].m_min[j];
        for(int j = 0; j < NUMDIMS; j++) query_file >> queries[i].m_max[j];
       
    }

    int found_count;
    int lin_count;
    clock_t q_time = 0;
    int rightestright = 0;
    int last_id = DATA_COUNT -1;
    clock_t insert_time, start_insert_time;
    clock_t start_query_time;
    for(int i = 0; i < query_size; i++){
        cerr << "QUERY " << i << "...\n";
        insert_time = 0.0;
        if(i % ratio == 0 && i != 0){
            for(int o =0; o < insert_count; o++){
                for(int j = 0; j < NUMDIMS; j++) data_file >> min[j];
                for(int j = 0; j < NUMDIMS; j++) data_file >> max[j];
                this_id = ++last_id;
                start_insert_time = clock();
               #if defined(GS) || defined(GSS) || defined(GSM) || defined(GSQ)
                    tree.Insert_gradual(min, max, this_id);
                #elif defined(CS) || defined(CSM)
                    tree.Insert_full(min, max, this_id);
                #else
                    tree.Insert(min, max, this_id);
                #endif
                insert_time += clock() - start_insert_time;
                data.push_back(Rect(min, max));
            }
        }
        start_query_time = clock();
         #if defined(GS) || defined(GSS) || defined(GSQ) || defined(GSM) || defined(CS) || defined(CSM)
            found_count = tree.QueryAdaptive(queries[i].m_min, queries[i].m_max);
            #ifdef stats
                int rightestright = tree.getRightestRight();
                int height = tree.TreeHeight();
                stats_file << "i " << rightestright << " " << count_total_holes << " " << height << " " << count_internal_nodes << " " << count_regular_leaves << " " << count_irregular_leaves << endl;
            #endif
        #else
            found_count = tree.QueryAdaptive(queries[i].m_min, queries[i].m_max);
            #ifdef stats
                int rightestright = tree.getRightestRight();
                int height = tree.TreeHeight();
                stats_file << "i " << rightestright << " " << count_total_holes << " " << height << " " << count_internal_nodes << " " << count_regular_leaves << " " << count_irregular_leaves << endl;
            #endif
        #endif
        q_time = clock() - start_query_time;
        
        times_file << found_count << " " << (double)insert_time/CLOCKS_PER_SEC << " " << (double)q_time/CLOCKS_PER_SEC <<  "\n";

        #ifdef linearcheck
            lin_count = 0;
            for(int ii =0; ii < data.size(); ii++){
                if(Overlap(&(queries[i]), &(data[ii]))) {
                    #ifdef debug
                        cerr << ii << endl;
                    #endif
                    lin_count++;
                }
            }

            if(lin_count != found_count){
                cout << "WRONG ANSWER, TRAGEDY, DISASTER, ALL THAT..." << endl;
                cout << "Query \n";
                for(int ii = 0; ii < NUMDIMS; ii++)
                {
                    cout << queries[i].m_min[ii] << " " << queries[i].m_max[ii] << ",";
                }
                cout << "\n";
                return 3;
            }
        #endif
	    cout << found_count << endl;
    }
    query_file.close();
    cout << "DONE!!\n";
    return 0;
}

