#ifndef GLIDE_CQRIPPLE_H
#define GLIDE_CQRIPPLE_H

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>

using namespace std;

#include <utility>
#include <queue>
#include <string>
#include <chrono>
#include <string.h>
#include <valarray>
#include "valarray"

#define ASSERT assert // RTree uses ASSERT( condition )
#ifndef Min
#define Min std::min
#endif //Min
#ifndef Max
#define Max std::max
#endif //Max

#ifdef timebreak
    clock_t scan_pending_time = 0.0;
    clock_t search_time = 0.0;
    clock_t add_ltas_time = 0.0;
    clock_t ripple_time = 0.0;

    // in search
    clock_t tba_assignment_time = 0.0;
    clock_t cracking_time = 0.0;
    clock_t traversedown_time = 0.0;

    // in ripple
    clock_t shuffle_time = 0.0;
#endif 

#define RTREE_TEMPLATE template<class DATATYPE, int TMAXNODES, int TMINNODES, int TMAXDATAPOINTS, int TMINDATAPOINTS>
#define RTREE_QUAL RTree<DATATYPE, TMAXNODES, TMINNODES, TMAXDATAPOINTS, TMINDATAPOINTS>

#define RTREE_DONT_USE_MEMPOOLS // This version does not contain a fixed memory allocator, fill in lines with EXAMPLE to implement one.

#ifdef data_size6M
    #define DATA_COUNT 6000000
#elif data_size3_8M
    #define DATA_COUNT 3800000
#elif data_size6_4M
    #define DATA_COUNT 6400000
#elif data_size6_8M
    #define DATA_COUNT 6800000
#elif data_size7_2M
    #define DATA_COUNT 7200000
#elif data_size7_6M
    #define DATA_COUNT 7600000
#elif data_size7_92M
    #define DATA_COUNT 79200000
#elif data_size9_5M
    #define DATA_COUNT 9500000
#elif data_size500K
    #define DATA_COUNT 500000
#elif data_size14_25M
    #define DATA_COUNT 14250000
#elif data_size15_2M
    #define DATA_COUNT 15200000
#elif data_size16_15M
    #define DATA_COUNT 16150000
#elif data_size17_1M
    #define DATA_COUNT 17100000
#elif data_size18_05M
    #define DATA_COUNT 18050000
#elif data_size18_81M
    #define DATA_COUNT 18810000
#elif data_size19M
    #define DATA_COUNT 19000000
#elif data_size12M
    #define DATA_COUNT 12000000
#elif data_size24M
    #define DATA_COUNT 24000000
#elif data_size48M
    #define DATA_COUNT 48000000
#elif data_size52_5M
    #define DATA_COUNT 52500000
#elif data_size66_5M
    #define DATA_COUNT 66500000
#elif data_size6_45M
    #define DATA_COUNT 6450000
#elif data_size8_17M
    #define DATA_COUNT 8170000
#elif data_size86_25M
    #define DATA_COUNT 86250000
#elif data_size109_25M
    #define DATA_COUNT 109250000
#endif


#define INITIAL_HOLES 0

#define NUMDIMS 2

static float m_data_arr_mins[DATA_COUNT * 3][NUMDIMS];
static float m_data_arr_maxes[DATA_COUNT * 3][NUMDIMS];
static int m_data_arr_ids[DATA_COUNT * 3];

static float m_pending_insertions_mins[DATA_COUNT][NUMDIMS];
static float m_pending_insertions_maxes[DATA_COUNT][NUMDIMS];
static int m_pending_insertions_ids[DATA_COUNT];

static float tba_mins[DATA_COUNT][NUMDIMS];
static float tba_maxes[DATA_COUNT][NUMDIMS];
static int tba_ids[DATA_COUNT];


int partition(int *mapping, int low, int high, int offset)
{
    // Choosing the pivot
    int folan;
    int pivot = mapping[high - offset];
 
    // Index of smaller element and indicates
    // the right position of pivot found so far
    int i = (low - 1);
 
    for (int j = low; j <= high - 1; j++) {
 
        // If current element is smaller than the pivot
        if (mapping[j - offset] < pivot) {
 
            // Increment index of smaller element
            i++;
            for(folan = 0; folan< NUMDIMS; folan++){
                swap(tba_mins[i][folan], tba_mins[j][folan]);
                swap(tba_maxes[i][folan], tba_maxes[j][folan]);
            }
            swap(tba_ids[i], tba_ids[j]);
            swap(mapping[i-offset],mapping[j-offset]);
        }
    }
    for(folan = 0; folan< NUMDIMS; folan++){
        swap(tba_mins[(i + 1)][folan], tba_mins[high][folan]);
        swap(tba_maxes[(i + 1)][folan], tba_maxes[high][folan]);
    }
    swap(tba_ids[i + 1], tba_ids[high]);
    swap(mapping[i + 1 - offset], mapping[high - offset]);
    return (i + 1);
}


void quickSort(int *mapping, int low, int high, int offset)
{
    if (low < high) {
 
        // pi is partitioning index, arr[p]
        // is now at right place
        int pi = partition(mapping, low, high, offset);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(mapping, low, pi - 1, offset);
        quickSort(mapping, pi + 1, high, offset);
    }
}



/// \class RTree
/// Implementation of RTree, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use RTree<Object*, float, 3> myTree;
///
/// This modified, templated C++ version by Greg Douglas at Auran (http://www.auran.com)
///
/// DATATYPE Referenced data, should be int, void*, obj* etc. no larger than sizeof<void*> and simple type
/// float Type of element such as int or float
/// NUMDIMS Number of dimensions such as 2 or 3
/// float Type of element that allows fractional and large values such as float or double, for use in volume calcs
///
/// NOTES: Inserting and removing data requires the knowledge of its constant Minimal Bounding Rectangle.
///        This version uses new/delete for nodes, I recommend using a fixed size allocator for efficiency.
///        Instead of using a callback function for returned results, I recommend and efficient pre-sized, grow-only memory
///        array similar to MFC CArray or STL Vector for returning search query result.
///
template<class DATATYPE,
        int TMAXNODES = 8, int TMINNODES = TMAXNODES / 2, int TMAXDATAPOINTS = 256, int TMINDATAPOINTS = TMAXDATAPOINTS / 2>
                class RTree {
                protected:

                    struct Node;  // Fwd decl.  Used by other internal structs and iterator

                public:

                    // These constant must be declared after Branch and before Node struct
                    // Stuck up here for MSVC 6 compiler.  NSVC .NET 2003 is much happier.
                    enum {
                        MAXNODES = TMAXNODES,                         ///< Max elements in node
                        MINNODES = TMINNODES,                         ///< Min elements in node
                        MAXDATAPOINTS = TMAXDATAPOINTS,
                        MINDATAPOINTS = TMINDATAPOINTS
                    };



                public:

                    // initialize rtree from data_file
                    // also do it with static memory
                    RTree(std::string data_file_name);

                    virtual ~RTree();

                    /// Insert entry
                    /// \param a_min Min of bounding rect
                    /// \param a_max Max of bounding rect
                    /// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
                    void Insert(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId);


                    int QueryAdaptive(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    /// Remove all entries from tree
                    void RemoveAll();

                    /// Count the data elements in this container.  This is slow as no internal counter is maintained.
                    int Count();

                    // stats about the tree
                    int CountDataPoints();
                    int SumDataPointIDs();
                    int CountNodes();
                    int TreeHeight();
                    std::pair<int, int> CountLeaves();
                    void printLeafLinkedList();
                    void printIDs();
                    int getPendingCount(){return m_pending_insertions_count;}
                    int getRightestRight(){return m_rightest_leaf->m_R;}
                    void printLeafBoundsOfData(int data_index);
                    bool printLeafBoundsOfDataRec(Node *a_node, int data_index);
                    void findAllLeavesWithDataInThem(int data_index);
                    void lookForEmptyNodes();

                protected:

                    /// Minimal bounding rectangle (n-dimensional)
                    struct Rect {
                        Rect(){}
                        Rect(float mins[NUMDIMS], float maxes[NUMDIMS]){
                            for(int i=0; i< NUMDIMS; i++) {
                                m_min[i] = mins[i];
                                m_max[i] = maxes[i];
                            }
                        }
                        float m_min[NUMDIMS];                      ///< Min dimensions of bounding box
                        float m_max[NUMDIMS];                      ///< Max dimensions of bounding box
                    };

                    /// May be data or may be another subtree
                    /// The parents level determines this.
                    /// If the parents level is 0, then this is data
                    struct Branch {
                        Rect m_rect;                                  ///< Bounds
                        Node *m_child;                                ///< Child node
                        DATATYPE m_data;                              ///< Data Id
                        Branch(){}
                        Branch(Node *the_node, Rect the_rect, DATATYPE the_data){
                            m_child = the_node;
                            m_data = the_data;
                            for(int i=0; i< NUMDIMS; i++) {
                                m_rect.m_min[i] = the_rect.m_min[i];
                                m_rect.m_max[i] = the_rect.m_max[i];
                            }
                        }
                        bool operator<(const Branch& a) const{
                            return m_child->m_L < a.m_child->m_L;
                        }
                    };

                    /// Node for each branch level
                    struct Node {
                        bool IsInternalNode() { return (m_level > 0); } // Not a leaf, but a internal node
                        bool IsLeaf() { return (m_level == 0); } // A leaf, contains data
                        // we always use the isReg function to check regularity
                        // to make sure the holes are always correct, let's just check here
                        // bool isRegular() { return m_regular; } // leaf is regular is m_regular is true, irregular ow
                        // we also only ask leaves
                        bool isRegular() { return ((m_R - m_L - m_holes) <= MAXDATAPOINTS); } // leaf is regular is m_regular is true, irregular ow

                        // for debug
                        std::string name;
                        int m_count;                                  ///< Count
                        int m_level;                                  ///< Leaf is zero, others positive
                        // SHOULD NOT USE THIS M_REGULAR THING.
                        // IT IS CORRECT I THINK
                        // BUT USE THE ISREGULAR FUNCTION!!
                        bool m_regular;                               ///< Some leaves can be irregular. Leaf is irregular if this is false, regular ow
                        //        std::vector<Branch> m_branch;                 ///< Branch
                        Branch m_branch[MAXNODES];                    ///< Branch
                        int m_L;                                        ///< left side of data_points interval in m_data_arr
                        int m_R;                                        ///< right side of data_points interval in m_data_arr
                        Node* m_parent;
                        int m_id;

                        // pointer to leaf on the right side
                        Node* m_right_sibling;
                        // pointer to leaf on the left side
                        Node* m_left_sibling;
                        // number of empty spaces at the front of this leaf's L.
                        int m_holes;
                        // MARCH 2023


                        Node(){}
                    };

                    /// Variables for finding a split partition
                    struct PartitionVars {
                        enum {
                            NOT_TAKEN = -1
                        }; // indicates that position

                        int m_partition[MAXDATAPOINTS + 1];
                        int m_total;
                        int m_minFill;
                        int m_count[2];
                        Rect m_cover[2];
                        float m_area[2];

                        Branch m_branchBuf[MAXDATAPOINTS + 1];
                        int m_branchCount;
                        Rect m_coverSplit;
                        float m_coverSplitArea;
                    };

                    typedef std::pair<bool,bool> bools;
                    typedef std::pair<std::vector<Branch>, std::vector<int>> tbd_vectors;
                    typedef std::pair<int, int> position;

                    struct LTA_v2{
                        Node* this_leaf;
                        int Ls[2*NUMDIMS + 1 + 1];
                        int Rs[2*NUMDIMS + 1 + 1];
                        int holes[2*NUMDIMS + 1 + 1];
                        float crack_covers_mins[2*NUMDIMS + 1 + 1][NUMDIMS];
                        float crack_covers_maxes[2*NUMDIMS + 1 + 1][NUMDIMS];
                        int how_many_created = 0;
                        int qp_L;
                        int tba_start;
                        int tba_count;

                    };

                    typedef vector<LTA_v2> LeavesToAdd_v2;

                    struct overlapping_leaf_tbas{
                        Node *this_leaf;
                        int this_tba_start;
                        int this_tba_count;
                        overlapping_leaf_tbas(Node* l, int tba_start, int tba_count){
                            this_leaf = l;
                            this_tba_start = tba_start;
                            this_tba_count = tba_count;
                        }
                    };


                    Node *AllocNode();

                    void FreeNode(Node *a_node);

                    void InitNode(Node *a_node);

                    void InitRect(Rect *a_rect);

                    bool InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level, Rect &a_coverRect, Rect &a_newCoverRect);

                    bool InsertRect(const Branch &a_branch, Node **a_root, int a_level);

                    void Insert_anylevel(const Branch &a_branch, Node *start, int a_level);
                    
                    int getNodeBranchIndex(Node* a_node);

                    Rect NodeCover(Node *a_node);

                    Rect LeafCover(int L, int R);

                    bool AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode);

                    bool AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect);

                    void DisconnectBranch(Node *a_node, int a_index);

                    int PickBranch(const Rect *a_rect, Node *a_node);

                    static Rect CombineRect(const Rect *a_rectA, const Rect *a_rectB);

                    void SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode);

                    void SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect);

                    float RectVolume(Rect *a_rect);

                    float CalcRectVolume(Rect *a_rect);

                    void GetBranches(Node *a_node, const Branch *a_branch, PartitionVars *a_parVars);

                    void ChoosePartition(PartitionVars *a_parVars, int a_minFill);

                    void LoadNodes(Node *a_nodeA, Node *a_nodeB, PartitionVars *a_parVars);

                    void InitParVars(PartitionVars *a_parVars, int a_maxRects, int a_minFill);

                    void PickSeeds(PartitionVars *a_parVars);

                    void Classify(int a_index, int a_group, PartitionVars *a_parVars);

                    bool Overlap(Rect *a_rectA, Rect *a_rectB) const;

                    bool Overlap(Rect *a_rectA, const float a_min[NUMDIMS], const float a_max[NUMDIMS]) const;

                    void TraverseDownTreeTilLeaf(Node* a_node, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves);

                    bool Search_2027(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves);

                    static bool compareByL(const Branch a, const Branch b);

                    void Add_ltas_v3(LeavesToAdd_v2 *all_lta, vector<overlapping_leaf_tbas> *overlapping_leaves);                    
                    
                    void deleteEmptyNodesUp(Node* start);
                    
                    bool is_leaf_in_OL(Node *a_leaf, vector<overlapping_leaf_tbas> *overlapping_leaves);

                    void shuffle_right(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k);
                    
                    void shuffle_left(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k);
                    
                    void ripple_v2(vector<overlapping_leaf_tbas> *overlapping_leaves);

                    void Reset();

                    void RemoveAllRec(Node *a_node);

                    void CopyBranch(Branch &current, const Branch &other);

                    // returns 2*axis + {0 if min or 1 if max}, for all choices ignore = -3
                    int ChooseAxisLargestSideMid(Rect node_cover, Rect query_cover);

                    int CrackOnAxisCompMin_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                   
                    int CrackOnAxisCompMax_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);

                    void swap(float &a, float &b);
                   
                    void swap_index(int i, int j);

                    int MyPartition_v45(int L, int R, float pivot_value, int axis, bool min_or_max, Rect *left_rect, Rect *right_rect);

                    bool CountDataPoints(Node *a_node, int &a_foundCount);
                    bool SumDataPointIDs(Node *a_node, int &a_foundCount);
                    bool CountNodes(Node *a_node, int &a_foundCount);
                    void PrintLeafSizesRec(Node *a_node, ofstream &myfile);
                    void lookForEmptyNodesRec(Node* a_node);
                    void findNodeCount(int node_id, Node* this_node);
                    
                    Node *m_root;                                    ///< Root of tree
                    
                    int last_id;
                    // the first leaf, the left most leaf
                    Node* m_leftest_leaf;
                    // last leaf, right most leaf
                    Node* m_rightest_leaf;

                    Rect root_cover;

                    int m_pending_insertions_count;
                };



RTREE_TEMPLATE
RTREE_QUAL::RTree(std::string data_file_name){
    std::ifstream data_file(data_file_name.c_str());
    if(data_file.is_open()){cout << data_file_name << endl;}

    data_file.clear();
    data_file.seekg(0, ios::beg);

    float min[NUMDIMS]; float max[NUMDIMS];

    for(int i = INITIAL_HOLES; i < DATA_COUNT+INITIAL_HOLES; i++){
        for(int j = 0; j < NUMDIMS; j++){
            data_file >> min[j];
        }
        for(int j = 0; j < NUMDIMS; j++){
            data_file >> max[j];
        }

        m_data_arr_ids[i] = i - INITIAL_HOLES;
        for (int axis = 0; axis < NUMDIMS; ++axis) {
            m_data_arr_mins[i][axis] = min[axis];
            m_data_arr_maxes[i][axis] = max[axis];
        }
    }

    data_file.close();

    m_root = AllocNode();
    m_root->m_level = 0;
    m_root->m_regular = false;
    m_root->m_parent = NULL;
    last_id = 0;
    m_root->m_id = last_id;

    last_id++;

    m_root->m_L = 0;
    m_root->m_R = DATA_COUNT + INITIAL_HOLES;
    m_root->m_holes = INITIAL_HOLES;
    m_root->m_count = m_root->m_R - m_root->m_L - m_root->m_holes;
    ASSERT(m_root->m_count == DATA_COUNT);

    root_cover = NodeCover(m_root);

    m_rightest_leaf = m_root;
    m_leftest_leaf = m_root;
    m_root->m_left_sibling = NULL;
    m_root->m_right_sibling = NULL;
    m_pending_insertions_count = 0;
}


RTREE_TEMPLATE
RTREE_QUAL::~RTree() {
    Reset(); // Free, or reset node memory
}


RTREE_TEMPLATE
void RTREE_QUAL::Insert(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId) {

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        m_pending_insertions_mins[m_pending_insertions_count][axis] = a_min[axis];
        m_pending_insertions_maxes[m_pending_insertions_count][axis] = a_max[axis];
    }

    m_pending_insertions_ids[m_pending_insertions_count] = a_dataId;

    m_pending_insertions_count ++;

}


RTREE_TEMPLATE
int RTREE_QUAL::CountDataPoints() {
    int foundCount = 0;
    CountDataPoints(m_root, foundCount);
    return foundCount;
}


RTREE_TEMPLATE
int RTREE_QUAL::CountNodes() {
    int foundCount = 0;
    CountNodes(m_root, foundCount);
    return foundCount;
}

RTREE_TEMPLATE
int RTREE_QUAL::TreeHeight() {
    return (m_root->m_level + 1);
}

RTREE_TEMPLATE
int RTREE_QUAL::SumDataPointIDs() {
    int foundCount = 0;
    SumDataPointIDs(m_root, foundCount);
    return foundCount;
}


RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive(const float *a_min, const float *a_max) {

    #ifdef timebreak
        scan_pending_time = 0.0;
        search_time = 0.0;
        add_ltas_time = 0.0;
        ripple_time = 0.0;

        // in search
        tba_assignment_time = 0.0;
        cracking_time = 0.0;
        traversedown_time = 0.0;

        // in ripple
        shuffle_time = 0.0;
    #endif 

    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    int tba_indexes[m_pending_insertions_count];
    int tba_count = 0;
    #ifdef timebreak
        auto start_time = clock();
    #endif
    for(int i = 0; i < m_pending_insertions_count; i++){
        if(Overlap(&rect, m_pending_insertions_mins[i], m_pending_insertions_maxes[i])){

            tba_indexes[tba_count] = i;
            tba_ids[tba_count] = m_pending_insertions_ids[i];

            for (int axis = 0; axis < NUMDIMS; ++axis) {
                tba_mins[tba_count][axis] = m_pending_insertions_mins[i][axis];
                tba_maxes[tba_count][axis] = m_pending_insertions_maxes[i][axis];
            }

            tba_count++;

            for (int axis = 0; axis < NUMDIMS; ++axis) {
                if(root_cover.m_min[axis] > m_pending_insertions_mins[i][axis])
                    {root_cover.m_min[axis] = m_pending_insertions_mins[i][axis];}
                if(root_cover.m_max[axis] < m_pending_insertions_maxes[i][axis])
                    {root_cover.m_max[axis] = m_pending_insertions_maxes[i][axis];}
            }
        }
    }

    #ifdef timebreak
        scan_pending_time = clock() - start_time;
    #endif

    for(int i = tba_count -1; i > -1; i--){
        if(tba_indexes[i] != m_pending_insertions_count - 1){
            for (int axis = 0; axis < NUMDIMS; ++axis) {
                m_pending_insertions_mins[tba_indexes[i]][axis] = m_pending_insertions_mins[m_pending_insertions_count -1][axis];
                m_pending_insertions_maxes[tba_indexes[i]][axis] = m_pending_insertions_maxes[m_pending_insertions_count - 1][axis];
            }
            m_pending_insertions_ids[tba_indexes[i]] = m_pending_insertions_ids[m_pending_insertions_count - 1];

            m_pending_insertions_count--;
        }
        else{
            m_pending_insertions_count--;
        }
    }

    int foundCount = 0;
    foundCount += tba_count;
    LeavesToAdd_v2 all_lta;

    vector<overlapping_leaf_tbas> ol;

    #ifdef timebreak
        start_time = clock();
    #endif
    Search_2027(m_root, &rect, foundCount, &all_lta, 0, tba_count, &ol);
    #ifdef timebreak
        search_time = clock() - start_time;
        start_time = clock();
    #endif
    Add_ltas_v3(&all_lta, &ol);
    #ifdef timebreak
        add_ltas_time = clock() - start_time;
        start_time = clock();
    #endif
    ripple_v2(&ol);
    #ifdef timebreak
        ripple_time = clock() - start_time;
        cout <<  (double)search_time/CLOCKS_PER_SEC << " " << (double)tba_assignment_time/CLOCKS_PER_SEC << " " << (double)cracking_time/CLOCKS_PER_SEC << " " << (double)traversedown_time/CLOCKS_PER_SEC << " ";
        cout << (double)add_ltas_time/CLOCKS_PER_SEC << " " << (double)ripple_time/CLOCKS_PER_SEC << " " << (double)shuffle_time/CLOCKS_PER_SEC << endl;
    #endif
    return foundCount;
}


RTREE_TEMPLATE
int RTREE_QUAL::Count() {
    int count = 0;
    CountRec(m_root, count);

    return count;
}


RTREE_TEMPLATE
void RTREE_QUAL::RemoveAll() {
    // Delete all existing nodes
    Reset();

    m_root = AllocNode();
    m_root->m_level = 0;
}


RTREE_TEMPLATE
void RTREE_QUAL::Reset() {
#ifdef RTREE_DONT_USE_MEMPOOLS
    // Delete all existing nodes
    RemoveAllRec(m_root);
#else // RTREE_DONT_USE_MEMPOOLS
    // Just reset memory pools.  We are not using complex types
    // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

RTREE_TEMPLATE
void RTREE_QUAL::RemoveAllRec(Node *a_node) {

    if (a_node->IsInternalNode()) // This is an internal node in the tree
    {
        for (int index = 0; index < a_node->m_count; ++index) {
            RemoveAllRec(a_node->m_branch[index].m_child);
        }
    }
    FreeNode(a_node);
}



RTREE_TEMPLATE
typename RTREE_QUAL::Node *RTREE_QUAL::AllocNode() {
    Node *newNode;
#ifdef RTREE_DONT_USE_MEMPOOLS
    newNode = new Node;
#else // RTREE_DONT_USE_MEMPOOLS
    // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
InitNode(newNode);
return newNode;
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeNode(Node *a_node) {

#ifdef RTREE_DONT_USE_MEMPOOLS
delete a_node;
#else // RTREE_DONT_USE_MEMPOOLS
// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::InitNode(Node *a_node) {
    a_node->m_count = 0;
    a_node->m_level = -1;
    a_node->m_regular = true;
    a_node->m_id = last_id;
    last_id++;
    a_node->m_right_sibling = NULL;
    a_node->m_left_sibling = NULL;
    a_node->m_holes = 0;
}


RTREE_TEMPLATE
void RTREE_QUAL::InitRect(Rect *a_rect) {
    for (int index = 0; index < NUMDIMS; ++index) {
        a_rect->m_min[index] = (float) 0;
        a_rect->m_max[index] = (float) 0;
    }
}


// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level, Rect &a_coverRect, Rect &a_newCoverRect) {
    // recurse until we reach the correct level for the new record. data records
    // will always be called with a_level == 0 (leaf)
    if (a_node->m_level > a_level) {
        // Still above level for insertion, go down tree recursively
        Node *otherNode;

        // find the optimal branch for this record
        int index = PickBranch(&a_branch.m_rect, a_node);

        // recursively insert this record into the picked branch
        bool childWasSplit = InsertRectRec(a_branch, a_node->m_branch[index].m_child, &otherNode, a_level, a_coverRect, a_newCoverRect);

        if (!childWasSplit) {
            // Child was not split. Merge the bounding box of the new record with the
            // existing bounding box
            a_node->m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));
            return false;
        } else {
            // Child was split. The old branches are now re-partitioned to two nodes
            // so we have to re-calculate the bounding boxes of each node
            a_node->m_branch[index].m_rect = Rect(a_coverRect.m_max, a_coverRect.m_max);
            Branch branch;
            branch.m_child = otherNode;
            branch.m_rect = Rect(a_newCoverRect.m_min, a_newCoverRect.m_max);

            // The old node is already a child of a_node. Now add the newly-created
            // node to a_node as well. a_node might be split because of that.
            return AddBranch(&branch, a_node, a_newNode, a_coverRect, a_newCoverRect);
        }
    }
    else if (a_node->m_level == a_level) {
        // We have reached level for insertion. Add rect, split if necessary
        return AddBranch(&a_branch, a_node, a_newNode, a_coverRect, a_newCoverRect);
    } else {
        // Should never occur
        return false;
    }
}


// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
//
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRect(const Branch &a_branch, Node **a_root, int a_level) {

    Node *newNode;

    Rect rect1, rect2;

    if (InsertRectRec(a_branch, *a_root, &newNode, a_level, rect1, rect2))  // Root split
        {
        // Grow tree taller and new root
        Node *newRoot = AllocNode();
        newRoot->m_level = (*a_root)->m_level + 1;
        newRoot->m_parent = NULL;

        Branch branch;

        // add old root node as a child of the new root
        branch.m_rect = Rect(rect1.m_min, rect1.m_max);
        branch.m_child = *a_root;
        AddBranch(&branch, newRoot, NULL);

        // add the split node as a child of the new root
        branch.m_rect = Rect(rect2.m_min, rect2.m_max);
        branch.m_child = newNode;
        AddBranch(&branch, newRoot, NULL);

        // set the new root as the root node
        *a_root = newRoot;

        return true;
        }

    return false;
}


// Find the smallest rectangle that includes all rectangles in branches of a node.
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::NodeCover(Node *a_node) {
    Rect rect;
    Rect alaki_rect;
    if( !a_node->IsLeaf()) {
        rect = a_node->m_branch[0].m_rect;
        for (int index = 1; index < a_node->m_count; ++index) {
            rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
        }
    } else{
        rect = Rect(m_data_arr_mins[a_node->m_L], m_data_arr_maxes[a_node->m_L]);
        for (int index = a_node->m_L + 1; index < a_node->m_R; ++index) {
            alaki_rect = Rect(m_data_arr_mins[index], m_data_arr_maxes[index]);
            rect = CombineRect(&rect, &(alaki_rect));
        }
    }
    return rect;
}


RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::LeafCover(int L, int R){
    Rect rect = Rect(m_data_arr_mins[L], m_data_arr_maxes[L]);
    for (int index = L + 1; index < R; ++index) {
        for(int i = 0; i < NUMDIMS; i++){
            if(rect.m_min[i] > m_data_arr_mins[index][i]) rect.m_min[i] = m_data_arr_mins[index][i];
            if(rect.m_max[i] < m_data_arr_maxes[index][i]) rect.m_max[i] = m_data_arr_maxes[index][i];
        }

    }
    return rect;
}


RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode) {
    if(a_node->IsInternalNode()){
        // regular internal node
        if (a_node->m_count < MAXNODES) // Split won't be necessary
        {
            a_node->m_branch[a_node->m_count] = *a_branch;
            (a_branch->m_child)->m_parent = a_node;
            ++a_node->m_count;
            return false;
        }
        else{
            // only called on internal nodes!:)
            SplitNode(a_node, a_branch, a_newNode);
            return true;
        }
    }
    else
    {
        // it's a leaf, so we have to check regularity
        if (a_node->isRegular()) 
        {
            cout << "INSERTING INTO REGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
        else
        {
            cout << "INSERTING INTO IRREGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
    }
}


RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect) {
    if(a_node->IsInternalNode()){
        // regular internal node
        if (a_node->m_count < MAXNODES) // Split won't be necessary
        {
            a_node->m_branch[a_node->m_count] = *a_branch;
            (a_branch->m_child)->m_parent = a_node;
            ++a_node->m_count;
            return false;
        }
        else{
            // only called on regular nodes!:)
            SplitNode(a_node, a_branch, a_newNode, a_coverRect, a_newCoverRect);
            return true;
        }
    }
    else
    {
        // it's a leaf, so we have to check regularity
        if (a_node->isRegular()) 
        {
            cout << "INSERTING INTO REGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
        else
        {
            cout << "INSERTING INTO IRREGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
    }
}


// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed
RTREE_TEMPLATE
void RTREE_QUAL::DisconnectBranch(Node *a_node, int a_index) {
    if(a_node->m_level == 1){
        if((a_node->m_branch[a_index]).m_child == m_leftest_leaf){
            // is head
            m_leftest_leaf = ((a_node->m_branch[a_index]).m_child)->m_right_sibling;
        }
        else{
            (((a_node->m_branch[a_index]).m_child)->m_left_sibling)->m_right_sibling = ((a_node->m_branch[a_index]).m_child)->m_right_sibling;
        }
        if((a_node->m_branch[a_index]).m_child == m_rightest_leaf){
            m_rightest_leaf = ((a_node->m_branch[a_index]).m_child)->m_left_sibling;
        }
        else{
            (((a_node->m_branch[a_index]).m_child)->m_right_sibling)->m_left_sibling = ((a_node->m_branch[a_index]).m_child)->m_left_sibling;
        }
    }
    // Remove element by swapping with the last element to prevent gaps in array
    a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];
    a_node->m_count -= 1;
    int max_children;
    if(a_node->IsLeaf()) { max_children = MAXDATAPOINTS;} else {max_children = MAXNODES;}
    if (a_node->m_count <= max_children) a_node->m_regular = true;
}


// Pick a branch.  Pick the one that will need the smallest increase
// in area to accomodate the new rectangle.  This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when searching.
RTREE_TEMPLATE
int RTREE_QUAL::PickBranch(const Rect *a_rect, Node *a_node) {

    bool firstTime = true;
    float increase;
    float bestIncr = (float) -1;
    float area;
    float bestArea;
    int best = 0;
    Rect tempRect;
    for (int index = 0; index < a_node->m_count; ++index) {
        
        Rect *curRect = &a_node->m_branch[index].m_rect;
        area = CalcRectVolume(curRect);
        tempRect = CombineRect(a_rect, curRect);
        increase = CalcRectVolume(&tempRect) - area;
        if ((increase < bestIncr) || firstTime) {
            best = index;
            bestArea = area;
            bestIncr = increase;
            firstTime = false;
        } else if ((increase == bestIncr) && (area < bestArea)) {
            best = index;
            bestArea = area;
            bestIncr = increase;
        }
    }
    return best;
}


// Combine two rectangles into larger one containing both
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::CombineRect(const Rect *a_rectA, const Rect *a_rectB) {
    Rect newRect;

    for (int index = 0; index < NUMDIMS; ++index) {
        newRect.m_min[index] = Min(a_rectA->m_min[index], a_rectB->m_min[index]);
        newRect.m_max[index] = Max(a_rectA->m_max[index], a_rectB->m_max[index]);
    }

    return newRect;
}



// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode) {

    // Could just use local here, but member or external is faster since it is reused
    PartitionVars localVars;
    PartitionVars *parVars = &localVars;

    // Load all the branches into a buffer, initialize old node
    GetBranches(a_node, a_branch, parVars);

    // Find partition
    int min_fill;
    if(a_node->IsLeaf()) {min_fill = MINDATAPOINTS;} else {min_fill = MINNODES;}
    ChoosePartition(parVars, min_fill);

    // Create a new node to hold (about) half of the branches
    *a_newNode = AllocNode();
    (*a_newNode)->m_level = a_node->m_level;
    (*a_newNode)->name = "from split node";

    // Put branches from buffer into 2 nodes according to the chosen partition
    a_node->m_count = 0;
    //    a_node->m_branch.clear();
    LoadNodes(a_node, *a_newNode, parVars);
}


RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect) {
    // Could just use local here, but member or external is faster since it is reused
    PartitionVars localVars;
    PartitionVars *parVars = &localVars;

    // Load all the branches into a buffer, initialize old node
    GetBranches(a_node, a_branch, parVars);

    // Find partition
    int min_fill;
    if(a_node->IsLeaf()) {min_fill = MINDATAPOINTS;} else {min_fill = MINNODES;}
    ChoosePartition(parVars, min_fill);

    // Create a new node to hold (about) half of the branches
    *a_newNode = AllocNode();
    (*a_newNode)->m_level = a_node->m_level;

    // Put branches from buffer into 2 nodes according to the chosen partition
    a_node->m_count = 0;
    LoadNodes(a_node, *a_newNode, parVars);
    a_coverRect = (parVars->m_cover[0]);
    a_newCoverRect = (parVars->m_cover[1]);
}


// Calculate the n-dimensional volume of a rectangle
RTREE_TEMPLATE
float RTREE_QUAL::RectVolume(Rect *a_rect) {

    float volume = (float) (a_rect->m_max[0] - a_rect->m_min[0]);
    float len;

    for (int index = 1; index < NUMDIMS; ++index) {
        len = a_rect->m_max[index] - a_rect->m_min[index];
        volume *= (len);
    }
    ASSERT(volume >= (float) 0);

    return volume;
}


// Use one of the methods to calculate retangle volume
RTREE_TEMPLATE
float RTREE_QUAL::CalcRectVolume(Rect *a_rect) {
#ifdef RTREE_USE_SPHERICAL_VOLUME
    return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
#else // RTREE_USE_SPHERICAL_VOLUME
return RectVolume(a_rect); // Faster but can cause poor merges
#endif // RTREE_USE_SPHERICAL_VOLUME
}


// Load branch buffer with branches from full node plus the extra branch.
RTREE_TEMPLATE
void RTREE_QUAL::GetBranches(Node *a_node, const Branch *a_branch, PartitionVars *a_parVars) {
    int max_children;
    if(a_node->IsLeaf()){max_children = MAXDATAPOINTS;} else {max_children = MAXNODES;}
    for (int index = 0; index < max_children; ++index) {
        if(a_node->IsLeaf()){
            a_parVars->m_branchBuf[index] = Branch(NULL, Rect(m_data_arr_mins[a_node->m_L + index], m_data_arr_maxes[a_node->m_L + index]), m_data_arr_ids[a_node->m_L + index]);
        } else {
            a_parVars->m_branchBuf[index] = a_node->m_branch[index];
        }
    }
    a_parVars->m_branchBuf[max_children] = *a_branch;
    a_parVars->m_branchCount = max_children + 1;

    a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
    for (int index = 1; index < max_children + 1; ++index) {
        a_parVars->m_coverSplit = CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
    }
    a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);
}


// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_QUAL::ChoosePartition(PartitionVars *a_parVars, int a_minFill) {

    float biggestDiff;
    int group, chosen = 0, betterGroup = 0;

    InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
    PickSeeds(a_parVars);

    while (((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
    && (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill))
    && (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill))) {
        biggestDiff = (float) -1;
        for (int index = 0; index < a_parVars->m_total; ++index) {
            if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index]) {
                Rect *curRect = &a_parVars->m_branchBuf[index].m_rect;
                Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
                Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);
                float growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];
                float growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];
                float diff = growth1 - growth0;
                if (diff >= 0) {
                    group = 0;
                } else {
                    group = 1;
                    diff = -diff;
                }

                if (diff > biggestDiff) {
                    biggestDiff = diff;
                    chosen = index;
                    betterGroup = group;
                } else if ((diff == biggestDiff) && (a_parVars->m_count[group] < a_parVars->m_count[betterGroup])) {
                    chosen = index;
                    betterGroup = group;
                }
            }
        }
        Classify(chosen, betterGroup, a_parVars);
    }

    // If one group too full, put remaining rects in the other
    if ((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total) {
        if (a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill) {
            group = 1;
        } else {
            group = 0;
        }
        for (int index = 0; index < a_parVars->m_total; ++index) {
            if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index]) {

                Classify(index, group, a_parVars);
            }
        }
    }
}


// Copy branches from the buffer into two nodes according to the partition.
RTREE_TEMPLATE
void RTREE_QUAL::LoadNodes(Node *a_nodeA, Node *a_nodeB, PartitionVars *a_parVars) {

    for (int index = 0; index < a_parVars->m_total; ++index) {

        int targetNodeIndex = a_parVars->m_partition[index];
        Node *targetNodes[] = {a_nodeA, a_nodeB};

        // It is assured that AddBranch here will not cause a node split.
        bool nodeWasSplit = AddBranch(&a_parVars->m_branchBuf[index], targetNodes[targetNodeIndex], NULL);
    }
}


// Initialize a PartitionVars structure.
RTREE_TEMPLATE
void RTREE_QUAL::InitParVars(PartitionVars *a_parVars, int a_maxRects, int a_minFill) {

    a_parVars->m_count[0] = a_parVars->m_count[1] = 0;
    a_parVars->m_area[0] = a_parVars->m_area[1] = (float) 0;
    a_parVars->m_total = a_maxRects;
    a_parVars->m_minFill = a_minFill;
    for (int index = 0; index < a_maxRects; ++index) {
        a_parVars->m_partition[index] = PartitionVars::NOT_TAKEN;
    }
}


RTREE_TEMPLATE
void RTREE_QUAL::PickSeeds(PartitionVars *a_parVars) {
    int seed0 = 0, seed1 = 0;
    float worst, waste;
    float area[a_parVars->m_total + 1];

    for (int index = 0; index < a_parVars->m_total; ++index) {
        area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
    }


    worst = -1 * a_parVars->m_coverSplitArea - 1;
    for (int indexA = 0; indexA < a_parVars->m_total - 1; ++indexA) {
        for (int indexB = indexA + 1; indexB < a_parVars->m_total; ++indexB) {
            Rect oneRect = CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect);
            waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
            if (waste > worst) {
                worst = waste;
                seed0 = indexA;
                seed1 = indexB;
            }
        }
    }
    Classify(seed0, 0, a_parVars);
    Classify(seed1, 1, a_parVars);
}


// Put a branch in one of the groups.
RTREE_TEMPLATE
void RTREE_QUAL::Classify(int a_index, int a_group, PartitionVars *a_parVars) {
    a_parVars->m_partition[a_index] = a_group;

    // Calculate combined rect
    if (a_parVars->m_count[a_group] == 0) {
        a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
    } else {
        a_parVars->m_cover[a_group] = CombineRect(&a_parVars->m_branchBuf[a_index].m_rect,
                                                  &a_parVars->m_cover[a_group]);
    }

    // Calculate volume of combined rect
    a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);

    ++a_parVars->m_count[a_group];
}



// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Rect *a_rectA, Rect *a_rectB) const {

    for (int index = 0; index < NUMDIMS; ++index) {
        if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
        a_rectB->m_min[index] > a_rectA->m_max[index]) {
            return false;
        }
    }
    return true;
}

// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Rect *a_rectA, const float a_min[NUMDIMS], const float a_max[NUMDIMS]) const {

    for (int index = 0; index < NUMDIMS; ++index) {
        if (a_rectA->m_min[index] > a_max[index] ||
        a_min[index] > a_rectA->m_max[index]) {
            return false;
        }
    }
    return true;
}

RTREE_TEMPLATE
void RTREE_QUAL::TraverseDownTreeTilLeaf(Node *a_node, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves){
    
    if(!a_node->IsInternalNode()){

        if(tba_count > 0){
            overlapping_leaves->push_back(overlapping_leaf_tbas(a_node, tba_start, tba_count));
            if(a_node->m_parent != NULL) {
                int a_node_branch_index = getNodeBranchIndex(a_node);
                Rect trash_rect;
                for(int i =tba_start; i < tba_start + tba_count; i++){
                    trash_rect = (Rect(tba_mins[i], tba_maxes[i]));
                    (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&trash_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
                }
            }
        }
        return;
    }
    int *tba_to_children_counts = (int *)calloc(a_node->m_count, sizeof(int));
    int *tba_map = (int *)calloc(tba_count,sizeof(int));

    int a_node_branch_index;
    if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);
    Rect this_rect;
    for(int i = tba_start; i < tba_start + tba_count; i++){
        this_rect = (Rect(tba_mins[i], tba_maxes[i]));
        tba_map[i-tba_start] = PickBranch(&this_rect, a_node);

        tba_to_children_counts[tba_map[i-tba_start]]++;

        if(a_node->m_parent != NULL) {
            (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&this_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
        }
    }

    quickSort(tba_map, tba_start,tba_start + tba_count-1, tba_start);
    
    int *tba_to_children_acc_counts = (int *)calloc(a_node->m_count, sizeof(int));
    tba_to_children_acc_counts[0] = 0;
    for(int i = 1; i < a_node->m_count ; i++)
        tba_to_children_acc_counts[i] = tba_to_children_acc_counts[i-1] + tba_to_children_counts[i-1];

    for(int i =0; i < a_node->m_count; i++){
        if(tba_to_children_counts[i] > 0){
            ASSERT(0 <= tba_start + tba_to_children_acc_counts[i] && tba_start + tba_to_children_acc_counts[i] < DATA_COUNT);
            TraverseDownTreeTilLeaf(a_node->m_branch[i].m_child, tba_start + tba_to_children_acc_counts[i], tba_to_children_counts[i], overlapping_leaves);
        }
    }

    free(tba_to_children_counts);
    free(tba_map);
    free(tba_to_children_acc_counts);
}


RTREE_TEMPLATE
bool RTREE_QUAL::Search_2027(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves) {

    if (a_node->IsInternalNode()) {
        #ifdef timebreak
            clock_t start_time = clock();
        #endif

        int *tba_to_children_counts = (int *)calloc(a_node->m_count, sizeof(int));
        int *tba_map = (int *)calloc(tba_count,sizeof(int));

        int a_node_branch_index;
        if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);
        // int branch_choice;
        Rect this_rect;
        for(int i = tba_start; i < tba_start + tba_count; i++){
            this_rect = (Rect(tba_mins[i], tba_maxes[i]));

            tba_map[i-tba_start] = PickBranch(&this_rect, a_node);

            tba_to_children_counts[tba_map[i-tba_start]]++;

            if(a_node->m_parent != NULL) {
                (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&this_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
            }
        }

        #ifdef timebreak
            tba_assignment_time += (clock() - start_time); 
        #endif

        quickSort(tba_map, tba_start, tba_start + tba_count - 1, tba_start);
        int *tba_to_children_acc_counts = (int *)calloc(a_node->m_count, sizeof(int));
        tba_to_children_acc_counts[0] = 0;
        for(int i = 1; i < a_node->m_count; i++)
            tba_to_children_acc_counts[i] = tba_to_children_acc_counts[i-1] + tba_to_children_counts[i-1];

        for (int index = 0; index < a_node->m_count; ++index) {

            if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                ASSERT(0 <= tba_start + tba_to_children_acc_counts[index] && tba_start + tba_to_children_acc_counts[index] < DATA_COUNT);

                if (!Search_2027(a_node->m_branch[index].m_child, a_rect, a_foundCount, all_lta, tba_start + tba_to_children_acc_counts[index], tba_to_children_counts[index], overlapping_leaves)) {
                    return false;
                }
            }
            else{
                if(tba_to_children_counts[index] > 0) {
                    #ifdef timebreak
                        start_time = clock();
                    #endif
                    ASSERT(0 <= tba_start + tba_to_children_acc_counts[index] && tba_start + tba_to_children_acc_counts[index] < DATA_COUNT);
                    TraverseDownTreeTilLeaf(a_node->m_branch[index].m_child, tba_start + tba_to_children_acc_counts[index], tba_to_children_counts[index], overlapping_leaves);
                    #ifdef timebreak
                    traversedown_time += (clock() - start_time);
                    #endif
                }
            }
        }

        free(tba_to_children_counts);
        free(tba_map);
        free(tba_to_children_acc_counts);
        
        
    }
    else {
        if(a_node->isRegular()) {
            int a_node_branch_index;
            if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);   
            for (int index = a_node->m_L + a_node->m_holes; index < a_node->m_R; index++) {
                if (Overlap(a_rect, m_data_arr_mins[index], m_data_arr_maxes[index])) {
                    a_foundCount++;
                }
            }
            if(tba_count > 0) {
                overlapping_leaves->push_back(overlapping_leaf_tbas(a_node, tba_start, tba_count));
                if(a_node->m_parent != NULL) {
                    Rect new_rect = (a_node->m_parent)->m_branch[a_node_branch_index].m_rect;
                    Rect trash_rect;
                    for(int ii = tba_start; ii < tba_start + tba_count; ii++){
                        trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                        new_rect = CombineRect(&trash_rect, &new_rect);
                    }
                    ((a_node->m_parent)->m_branch[a_node_branch_index]).m_rect = new_rect;
                }
            }
        }
        else{

            #ifdef timebreak
                auto start_time = clock();
            #endif

            Rect node_cover;
            if(a_node->m_parent != NULL) {
                node_cover = ((a_node->m_parent)->m_branch[getNodeBranchIndex(a_node)]).m_rect;
            } else{
                node_cover = root_cover;
            }


            Rect this_piece_rect = Rect(node_cover.m_min, node_cover.m_max);
            int this_piece_L = a_node->m_L;
            int this_piece_R = a_node->m_R;

            int this_piece_holes = a_node->m_holes;


            int choice, mm, chosen_axis;

            LTA_v2 this_lta;
            this_lta.how_many_created =0;
            this_lta.this_leaf = a_node;
            Rect potential_query_piece = Rect(node_cover.m_min, node_cover.m_max); 
            Rect other_piece;


            int crack_index, filan, ashghal;
            int sum_of_choices = 0; int sum_of_mms = 0;

            int largest_piece_index = 0;
            int stochastic_crack_axis = 0; float stochastic_crack_pivot;
            float X, X2, X3;
            int stochastic_crack_index;

            for(int i =0; i < 2*NUMDIMS; i++){
                for(filan=0; filan < NUMDIMS; filan++) {
                    other_piece.m_min[filan] = potential_query_piece.m_min[filan];
                    other_piece.m_max[filan] = potential_query_piece.m_max[filan];
                }
                if(i == ((2 * NUMDIMS) - 1)){
                    choice = (((2*NUMDIMS) - 1) * NUMDIMS) - sum_of_choices;
                }
                else {
                    choice = ChooseAxisLargestSideMid(this_piece_rect, *a_rect);
                }
                if(choice < 0) break;
                if(sum_of_mms >= NUMDIMS) mm = 1;
                else if (sum_of_mms <= (-1*NUMDIMS)) mm = 0;
                else {
                    mm = choice % 2;
                    if(mm == 0) sum_of_mms++;
                    else sum_of_mms--;
                }

                chosen_axis = choice / 2;
                if(chosen_axis >= NUMDIMS || chosen_axis < 0) break;

                sum_of_choices += choice;

                other_piece.m_min[chosen_axis] = 21474836;
                other_piece.m_max[chosen_axis] = -21474836;

                if(mm == 0) {
                    potential_query_piece.m_min[chosen_axis] = 21474836;
                    crack_index = CrackOnAxisCompMax_v45( this_piece_L + this_piece_holes, this_piece_R, a_rect->m_min[chosen_axis],
                                                        chosen_axis, &other_piece, &potential_query_piece);
                }
                else{
                    potential_query_piece.m_max[chosen_axis] = -21474836;
                    crack_index = CrackOnAxisCompMin_v45( this_piece_L + this_piece_holes, this_piece_R, a_rect->m_max[chosen_axis],
                                                        chosen_axis, &potential_query_piece, &other_piece);
                }


                if((mm == 0 && crack_index != this_piece_L + this_piece_holes) ||
                (mm == 1 && this_piece_R != crack_index)){

                    // distribute holes

                    if(mm == 0){
                        
                        int half_the_holes = (int)(this_piece_holes / 2);
                        this_lta.Ls[this_lta.how_many_created] = this_piece_L;

                        if(this_piece_R != crack_index){
                            if(crack_index - this_piece_holes - this_piece_L >= half_the_holes)
                            {
                                this_lta.holes[this_lta.how_many_created] = this_piece_holes - half_the_holes;

                                // moving h/2 from the end of this piece to the start, to leave the rest of the holes to be distributed again
                                memcpy(m_data_arr_mins[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_mins[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_maxes[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + this_piece_holes - half_the_holes, m_data_arr_ids + crack_index - half_the_holes, half_the_holes*sizeof(int));

                                this_lta.Rs[this_lta.how_many_created] = crack_index - half_the_holes;

                                this_piece_L = crack_index - half_the_holes;
                                this_piece_holes = half_the_holes;
                            }
                            else
                            {
                                this_lta.holes[this_lta.how_many_created] = this_piece_holes - half_the_holes;
                                int how_many_to_move = crack_index - this_piece_holes - this_piece_L;

                                memcpy(m_data_arr_mins[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_mins[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_maxes[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + this_piece_holes - half_the_holes, m_data_arr_ids + crack_index - how_many_to_move, how_many_to_move*sizeof(int));

                                this_lta.Rs[this_lta.how_many_created] = crack_index - half_the_holes;
                                this_piece_L = crack_index - half_the_holes;
                                this_piece_holes = half_the_holes;

                            }

                        }else {
                            this_lta.Rs[this_lta.how_many_created] = crack_index;
                            this_lta.holes[this_lta.how_many_created] = this_piece_holes;

                            this_piece_holes = 0; 
                            this_piece_L = crack_index;
                        }

                    }
                    else{
                        int half_the_holes = (int)(this_piece_holes / 2);
                        this_lta.Rs[this_lta.how_many_created] = this_piece_R;

                        if(crack_index != this_piece_L + this_piece_holes){
                            if(crack_index - this_piece_holes - this_piece_L >= half_the_holes)
                            {
                                this_lta.holes[this_lta.how_many_created] = half_the_holes;

                                // moving h/2 from the end of this piece to the start, to leave the rest of the holes to be distributed again
                                memcpy(m_data_arr_mins[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_mins[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_maxes[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + this_piece_holes - half_the_holes, m_data_arr_ids + crack_index - half_the_holes, half_the_holes*sizeof(int));

                                this_lta.Ls[this_lta.how_many_created] = crack_index - half_the_holes;

                                this_piece_R = crack_index - half_the_holes;
                                this_piece_holes -= half_the_holes;
                            }
                            else
                            {
                                this_lta.holes[this_lta.how_many_created] = this_piece_holes - half_the_holes;
                                int how_many_to_move = crack_index - this_piece_holes - this_piece_L;

                                // moving h/2 from the end of this piece to the start, to leave the rest of the holes to be distributed again
                                memcpy(m_data_arr_mins[this_piece_L + half_the_holes], m_data_arr_mins[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + half_the_holes], m_data_arr_maxes[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + half_the_holes, m_data_arr_ids + crack_index - how_many_to_move, how_many_to_move*sizeof(int));

                                this_lta.Ls[this_lta.how_many_created] = this_piece_L + half_the_holes + how_many_to_move;

                                this_piece_R = this_piece_L + half_the_holes + how_many_to_move;
                                this_piece_holes = half_the_holes;
                            }
                        } else {
                            ASSERT(crack_index - this_piece_L == this_piece_holes);
                            this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                            this_lta.holes[this_lta.how_many_created] = this_piece_holes;

                            this_piece_R = crack_index;
                        }
                    }
                    this_piece_rect = Rect(potential_query_piece.m_min, potential_query_piece.m_max);

                    for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                        this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = other_piece.m_min[ashghal];
                        this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = other_piece.m_max[ashghal];
                    }

                    this_lta.how_many_created++;
                }
                else{
                    if( mm == 0 ) {this_piece_rect.m_min[chosen_axis] = a_rect->m_min[chosen_axis]; }
                    else {this_piece_rect.m_max[chosen_axis] = a_rect->m_max[chosen_axis];}
                }

                if((this_piece_R - this_piece_L - this_piece_holes) <= MAXDATAPOINTS) {
                    break;
                }

            }

            // scan last piece
            for(int folan = this_piece_L + this_piece_holes; folan < this_piece_R; folan++){
                if (Overlap(a_rect, m_data_arr_mins[folan], m_data_arr_maxes[folan])) {
                    a_foundCount++;
                }
            }
            this_lta.qp_L = this_piece_L;
            this_lta.tba_start = tba_start;
            this_lta.tba_count = tba_count;

            if(this_piece_L + this_piece_holes != this_piece_R) {
                this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                this_lta.Rs[this_lta.how_many_created] = this_piece_R;

                this_lta.holes[this_lta.how_many_created] = this_piece_holes;

                for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                    this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = potential_query_piece.m_min[ashghal];
                    this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = potential_query_piece.m_max[ashghal];
                }
                this_lta.how_many_created++;
            }
            else{
                this_lta.qp_L = this_lta.Ls[this_lta.how_many_created - 1];
            }
            for(ashghal = 1; ashghal < this_lta.how_many_created; ashghal++){
                if((this_lta.Rs[ashghal] - this_lta.Ls[ashghal]) > (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])) largest_piece_index = ashghal;
            }

            if(this_lta.Ls[largest_piece_index] != this_lta.qp_L){
                
                stochastic_crack_axis = 0;
                for(ashghal = 1; ashghal < NUMDIMS; ashghal++){
                    if((this_lta.crack_covers_maxes[largest_piece_index][ashghal] - this_lta.crack_covers_mins[largest_piece_index][ashghal])  > (this_lta.crack_covers_maxes[largest_piece_index][stochastic_crack_axis] - this_lta.crack_covers_mins[largest_piece_index][stochastic_crack_axis]) ){
                        stochastic_crack_axis = ashghal;
                    }
                }

                X = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position
                X2 = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position
                X3 = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position

                // make X the median of the samples
                if (X2 > X) swap(X, X2);
                if (X > X3) swap(X, X3);
                if (X2 > X) swap(X, X2);
                stochastic_crack_pivot = X;

                Rect l_rect = Rect(this_lta.crack_covers_mins[largest_piece_index], this_lta.crack_covers_maxes[largest_piece_index]);
                Rect r_rect = Rect(this_lta.crack_covers_mins[largest_piece_index], this_lta.crack_covers_maxes[largest_piece_index]);
                l_rect.m_min[stochastic_crack_axis] = 21474836;
                l_rect.m_max[stochastic_crack_axis] = -21474836;
                r_rect.m_min[stochastic_crack_axis] = 21474836;

                stochastic_crack_index = CrackOnAxisCompMax_v45( this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index], this_lta.Rs[largest_piece_index], stochastic_crack_pivot,
                                                                stochastic_crack_axis, &l_rect, &r_rect);
                if(stochastic_crack_index != this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index] && stochastic_crack_index != this_lta.Rs[largest_piece_index]){
                    // distribute the holes
                    int half_the_holes = (int)(this_lta.holes[largest_piece_index] / 2);

                    if(stochastic_crack_index - this_lta.holes[largest_piece_index] - this_lta.Ls[largest_piece_index] >= half_the_holes)
                    {
                        this_lta.Ls[this_lta.how_many_created] = stochastic_crack_index - half_the_holes;
                        this_lta.Rs[this_lta.how_many_created] = this_lta.Rs[largest_piece_index];
                        this_lta.holes[this_lta.how_many_created] = half_the_holes;
                        for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                            this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = r_rect.m_min[ashghal];
                            this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = r_rect.m_max[ashghal];
                        }

                        this_lta.how_many_created++;
                        
                        this_lta.holes[largest_piece_index] = this_lta.holes[largest_piece_index] - half_the_holes;

                        memcpy(m_data_arr_mins[this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index]], m_data_arr_mins[stochastic_crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                        memcpy(m_data_arr_maxes[this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index]], m_data_arr_maxes[stochastic_crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                        memcpy(m_data_arr_ids + this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index], m_data_arr_ids + stochastic_crack_index - half_the_holes, half_the_holes*sizeof(int));

                        this_lta.Rs[largest_piece_index] = stochastic_crack_index - half_the_holes;

                        for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                            this_lta.crack_covers_mins[largest_piece_index][ashghal] = l_rect.m_min[ashghal];
                            this_lta.crack_covers_maxes[largest_piece_index][ashghal] = l_rect.m_max[ashghal];
                        }
                    }
                    else
                    {
                        this_lta.Ls[this_lta.how_many_created] = stochastic_crack_index - half_the_holes;
                        this_lta.Rs[this_lta.how_many_created] = this_lta.Rs[largest_piece_index];
                        this_lta.holes[this_lta.how_many_created] = half_the_holes;
                        for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                            this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = r_rect.m_min[ashghal];
                            this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = r_rect.m_max[ashghal];
                        }

                        this_lta.how_many_created++;

                        int old_holes = this_lta.holes[largest_piece_index];

                        int how_many_to_move = stochastic_crack_index - this_lta.holes[largest_piece_index] - this_lta.Ls[largest_piece_index];
                        this_lta.holes[largest_piece_index] = this_lta.holes[largest_piece_index] - half_the_holes;

                        memcpy(m_data_arr_mins[this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index]], m_data_arr_mins[stochastic_crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                        memcpy(m_data_arr_maxes[this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index]], m_data_arr_maxes[stochastic_crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                        memcpy(m_data_arr_ids + this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index], m_data_arr_ids + stochastic_crack_index - how_many_to_move, how_many_to_move*sizeof(int));

                        this_lta.Rs[largest_piece_index] = stochastic_crack_index - half_the_holes;
                    }
                }

           }
            // END OF STOCH 
            all_lta->push_back(this_lta);
            #ifdef timebreak
                cracking_time += (clock() - start_time);
            #endif
        }
    }
    return true; // Continue searching
}


RTREE_TEMPLATE
bool RTREE_QUAL::compareByL(const Branch a, const Branch b){
    if(a.m_child->m_L < b.m_child->m_L) return true;
    return false;
}

RTREE_TEMPLATE
void RTREE_QUAL::Add_ltas_v3(LeavesToAdd_v2 *all_lta, vector<overlapping_leaf_tbas> *overlapping_leaves){

    Node* this_leaf_older_brother; Node* this_leaf_younger_brother;

    Rect trash_rect;

    for(int i = 0; i < all_lta->size(); i++){
        if((all_lta->at(i)).how_many_created == 1) 
        {
            if(all_lta->at(i).tba_count > 0){
                overlapping_leaves->push_back(overlapping_leaf_tbas(all_lta->at(i).this_leaf, all_lta->at(i).tba_start, all_lta->at(i).tba_count));
                Rect new_rect; int branch_index;
                Rect trash_rect;
                if((all_lta->at(i).this_leaf)->m_parent != NULL){
                    branch_index = getNodeBranchIndex(all_lta->at(i).this_leaf);
                    new_rect = ((all_lta->at(i).this_leaf)->m_parent)->m_branch[branch_index].m_rect;
                
                    for(int ii = (all_lta->at(i).tba_start); ii < (all_lta->at(i).tba_start) + (all_lta->at(i).tba_count); ii++){
                            trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                            new_rect = CombineRect(&(new_rect), &trash_rect);
                        }

                    ((all_lta->at(i).this_leaf)->m_parent)->m_branch[branch_index].m_rect = new_rect;
                }
            }

            continue;
        }

        vector<Branch> tba((all_lta->at(i)).how_many_created);

        int old_id = (all_lta->at(i).this_leaf)->m_id;
        
        if((all_lta->at(i).this_leaf)->m_parent != NULL){
            Node* parent_of_start = (all_lta->at(i).this_leaf)->m_parent;
            int branch_index = getNodeBranchIndex((all_lta->at(i).this_leaf));
            this_leaf_older_brother = (all_lta->at(i).this_leaf)->m_left_sibling;
            this_leaf_younger_brother = (all_lta->at(i).this_leaf)->m_right_sibling;
            DisconnectBranch(parent_of_start, branch_index);
            FreeNode(all_lta->at(i).this_leaf);
            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                Branch this_branch;
                this_branch.m_child = AllocNode();
                this_branch.m_child->m_L = (all_lta->at(i)).Ls[j];
                this_branch.m_child->m_R = (all_lta->at(i)).Rs[j];

                this_branch.m_child->m_holes = (all_lta->at(i)).holes[j];

                this_branch.m_child->m_level = 0;
                this_branch.m_child->m_count = this_branch.m_child->m_R - this_branch.m_child->m_L - this_branch.m_child->m_holes;

                if(this_branch.m_child->m_count > MAXDATAPOINTS) {
                    this_branch.m_child->m_regular = false;
                    this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);
                }
                else {
                    this_branch.m_child->m_regular = true;
                    this_branch.m_rect = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                }
                
                // let's also extend the MBR here
                if(all_lta->at(i).qp_L == this_branch.m_child->m_L && (all_lta->at(i).tba_count) > 0){
                    for(int ii = (all_lta->at(i).tba_start); ii < (all_lta->at(i).tba_start) + (all_lta->at(i).tba_count); ii++){
                        trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                        this_branch.m_rect = CombineRect(&(this_branch.m_rect), &trash_rect);
                    }
                }
                tba[j] = this_branch;
            }

            std::sort(tba.begin(), tba.end());

            int tba_size;
            tba_size = (all_lta->at(i).tba_count);
            int hmc = ((all_lta->at(i)).how_many_created);
            tba[0].m_child->m_left_sibling = this_leaf_older_brother;
            if(this_leaf_older_brother != NULL) this_leaf_older_brother->m_right_sibling = tba[0].m_child;
            else m_leftest_leaf = tba.at(0).m_child;
            tba[hmc - 1].m_child->m_right_sibling = this_leaf_younger_brother;
            if(this_leaf_younger_brother != NULL) this_leaf_younger_brother->m_left_sibling = tba[hmc- 1].m_child;
            else m_rightest_leaf = tba[hmc-1].m_child;


            for(int j = 0; j < hmc; j++){
                if(j ==0){
                    tba.at(j).m_child->m_right_sibling = tba.at(j + 1).m_child;
                }
                else if(j == hmc - 1){
                    tba.at(j).m_child->m_left_sibling = tba.at(j-1).m_child;
                }
                else{
                    tba.at(j).m_child->m_right_sibling = tba.at(j + 1).m_child;
                    tba.at(j).m_child->m_left_sibling = tba.at(j-1).m_child;
                }

                ASSERT(tba.at(j).m_child->m_count > 0);

                Insert_anylevel(tba.at(j), parent_of_start, 1);

                if((all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L)){
                    if(all_lta->at(i).tba_count > 0)
                        overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).tba_start, all_lta->at(i).tba_count));
                }
            }
        }
        else{
            m_root->m_count = 0;
            m_root->m_level++;
            m_root->m_regular = true;
            m_root->m_parent = NULL;

            this_leaf_older_brother = NULL; this_leaf_younger_brother = NULL;

            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                Branch this_branch;
                this_branch.m_child = AllocNode();
                this_branch.m_child->m_L = (all_lta->at(i)).Ls[j];
                this_branch.m_child->m_R = (all_lta->at(i)).Rs[j];
                this_branch.m_child->m_holes = (all_lta->at(i)).holes[j];

                this_branch.m_child->m_level = 0;
                this_branch.m_child->m_count = this_branch.m_child->m_R - this_branch.m_child->m_L  - this_branch.m_child->m_holes;

                if(this_branch.m_child->m_count > MAXDATAPOINTS) {
                    this_branch.m_child->m_regular = false;
                    this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);
                }
                else {
                    this_branch.m_child->m_regular = true;
                    this_branch.m_rect = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                }

                // let's also extend the MBR here
                if(all_lta->at(i).qp_L == this_branch.m_child->m_L && (all_lta->at(i).tba_count) > 0){
                    Rect trash_rect;
                    for(int ii = (all_lta->at(i).tba_start); ii < (all_lta->at(i).tba_start) + (all_lta->at(i).tba_count); ii++){
                        trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                        this_branch.m_rect = CombineRect(&(this_branch.m_rect), &trash_rect);
                    }
                }
                tba[j] = this_branch;
            }

            std::sort(tba.begin(), tba.end());

            int hmc = ((all_lta->at(i)).how_many_created);
            tba[0].m_child->m_left_sibling = this_leaf_older_brother;
            if(this_leaf_older_brother != NULL) this_leaf_older_brother->m_right_sibling = tba[0].m_child;
            else m_leftest_leaf = tba.at(0).m_child;
            tba[hmc - 1].m_child->m_right_sibling = this_leaf_younger_brother;
            if(this_leaf_younger_brother != NULL) this_leaf_younger_brother->m_left_sibling = tba[hmc- 1].m_child;
            else m_rightest_leaf = tba[hmc-1].m_child;

            for(int j = 0; j < hmc; j++){
                if(j ==0){
                    tba.at(j).m_child->m_right_sibling = tba.at(j + 1).m_child;
                }
                else if(j == hmc - 1){
                    tba.at(j).m_child->m_left_sibling = tba.at(j-1).m_child;
                }
                else{
                    tba.at(j).m_child->m_right_sibling = tba.at(j + 1).m_child;
                    tba.at(j).m_child->m_left_sibling = tba.at(j-1).m_child;
                }
                ASSERT(tba.at(j).m_child->m_count > 0);
                InsertRect(tba.at(j), &m_root, 1);

                if( (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L)){
                    if(all_lta->at(i).tba_count > 0)
                        overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).tba_start, all_lta->at(i).tba_count));
                }
            }
            ASSERT(m_root->m_count == hmc);
        }
    }
}

RTREE_TEMPLATE
void RTREE_QUAL::deleteEmptyNodesUp(Node *start){
    if(start == NULL){
        return;
    }
    
    if(start->m_parent == NULL){
        // is the root. we don't do anything.
        return;
    }
    if(start->m_count > 0){
        return;
    }
    
    // it is empty
    Node* parent_of_start = start->m_parent;
    int branch_index = getNodeBranchIndex(start);
    DisconnectBranch(parent_of_start, branch_index);
    FreeNode(start);
    deleteEmptyNodesUp(parent_of_start);
}


RTREE_TEMPLATE
bool RTREE_QUAL::is_leaf_in_OL(Node *a_leaf, vector<overlapping_leaf_tbas> *overlapping_leaves){
    for(int folan = 0; folan < overlapping_leaves->size(); folan++){
        if(a_leaf->m_id == (overlapping_leaves->at(folan)).this_leaf->m_id ) return true;
    }
    return false;
}


RTREE_TEMPLATE
void RTREE_QUAL::shuffle_right(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k){

    if(k <= 0) {
        return;
    }

    int folan, filan;

    if(k <= this_leaf->m_holes){
        memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));

        this_leaf->m_L += k;
        (this_leaf->m_left_sibling)->m_R += k;
        this_leaf->m_holes -= k;
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        return;
    }

    else if(this_leaf->m_id == m_rightest_leaf->m_id){
        if(k > (this_leaf->m_R - this_leaf->m_L)){
            memcpy(m_data_arr_mins[this_leaf->m_L + k], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L + k], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L + k, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));
        }
        else{
            // % holes < k < size
            memcpy(m_data_arr_mins[this_leaf->m_R], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));            
            memcpy(m_data_arr_maxes[this_leaf->m_R], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));            
            memcpy(m_data_arr_ids + this_leaf->m_R, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (k - this_leaf->m_holes)*sizeof(int));

        }
        
        memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));

        this_leaf->m_L += k;
        this_leaf->m_R += (k - this_leaf->m_holes);
        (this_leaf->m_left_sibling)->m_R += k;
        this_leaf->m_holes = 0;
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        return;
    }
    // if this leaf is not in overlapping leaves
    else if(!is_leaf_in_OL(this_leaf, overlapping_leaves)){

        if(k >= (this_leaf->m_R - this_leaf->m_L)){
            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));
            m_pending_insertions_count += (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);

            memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, (this_leaf->m_R - this_leaf->m_L)*sizeof(int));

            (this_leaf->m_left_sibling)->m_R = this_leaf->m_R;

            Node *next_leaf = this_leaf->m_right_sibling; 
            int this_leaf_size = this_leaf->m_R - this_leaf->m_L;

            // remove this leaf from parent
            int branch_index = getNodeBranchIndex(this_leaf);
            Node* parent = this_leaf->m_parent;
            DisconnectBranch(parent, branch_index);
            FreeNode(this_leaf);
            deleteEmptyNodesUp(parent);
            shuffle_right(overlapping_leaves, next_leaf, this_tba_start + this_leaf_size, k - this_leaf_size);
        }
        else{
            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (k - this_leaf->m_holes)*sizeof(int));
            
            m_pending_insertions_count += (k - this_leaf->m_holes);
            
            memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));
            
            this_leaf->m_L += k;
            (this_leaf->m_left_sibling)->m_R += k;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            return;
        }
    
    }
    else{
        if(k >= (this_leaf->m_R - this_leaf->m_L)){
            Node * next_leaf = this_leaf->m_right_sibling;
            (this_leaf->m_left_sibling)->m_R += (this_leaf->m_R - this_leaf->m_L);
            int this_old_L = this_leaf->m_L, this_old_R = this_leaf->m_R;
            int last_filled_spot_in_arr = m_rightest_leaf->m_R;

            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));

            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
            this_leaf->m_L = last_filled_spot_in_arr;
            this_leaf->m_R = last_filled_spot_in_arr + folan;

            (this_leaf->m_left_sibling)->m_right_sibling = this_leaf->m_right_sibling;
            (this_leaf->m_right_sibling)->m_left_sibling = this_leaf->m_left_sibling;

            this_leaf->m_left_sibling = m_rightest_leaf;
            m_rightest_leaf->m_right_sibling = this_leaf;
            this_leaf->m_right_sibling = NULL;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            m_rightest_leaf = this_leaf;

            memcpy(m_data_arr_mins[this_old_L], tba_mins[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_old_L], tba_maxes[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_old_L, tba_ids + this_tba_start, (this_old_R - this_old_L)*sizeof(int));
            shuffle_right(overlapping_leaves, next_leaf, this_tba_start + this_old_R - this_old_L, k - (this_old_R - this_old_L));
        }
        else{
        float* temp_mins;
        temp_mins = (float*)malloc((k - this_leaf->m_holes) * NUMDIMS * sizeof(float));

        float* temp_maxes;
        temp_maxes = (float*)malloc((k - this_leaf->m_holes)  * NUMDIMS* sizeof(float));

        int* temp_ids = (int*)malloc((k - this_leaf->m_holes) * sizeof(int));

        memcpy(temp_mins,  m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(temp_maxes,  m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(temp_ids, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (k - this_leaf->m_holes)*sizeof(int));
        
        memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));

        memcpy(tba_mins[this_tba_start], temp_mins, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(tba_maxes[this_tba_start], temp_maxes, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(tba_ids + this_tba_start, temp_ids, (k - this_leaf->m_holes)*sizeof(int));
        
        this_leaf->m_L += k;
        (this_leaf->m_left_sibling)->m_R += k;
        k -= this_leaf->m_holes;
        this_leaf->m_holes = 0;
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        shuffle_right(overlapping_leaves, this_leaf->m_right_sibling, this_tba_start, k);

        free(temp_ids);
        free(temp_mins);
        free(temp_maxes);
        }
    }

}


RTREE_TEMPLATE
void RTREE_QUAL::shuffle_left(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k){
    if(k <= 0) return;

    int folan, filan;
    if(k <= this_leaf->m_holes){
        int src = max(this_leaf->m_R - k, this_leaf->m_L + this_leaf->m_holes);
        int how_many_to_move = min(k, this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
        int dest = min(this_leaf->m_L + this_leaf->m_holes - how_many_to_move, this_leaf->m_R - k - how_many_to_move);

        memcpy(m_data_arr_mins[dest], m_data_arr_mins[src], NUMDIMS*how_many_to_move*sizeof(float));
        memcpy(m_data_arr_maxes[dest], m_data_arr_maxes[src], NUMDIMS*how_many_to_move*sizeof(float));
        memcpy(m_data_arr_ids + dest, m_data_arr_ids + src, how_many_to_move*sizeof(int));

        memcpy(m_data_arr_mins[this_leaf->m_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_R - k, tba_ids + this_tba_start, k*sizeof(int));
        this_leaf->m_holes = dest - this_leaf->m_L;
        this_leaf->m_R -= k;
        (this_leaf->m_right_sibling)->m_L -= k;
        
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        (this_leaf->m_right_sibling)->m_holes = 0;
        return;
    }
    
    else if(this_leaf->m_id == m_leftest_leaf->m_id){
        if(k > (this_leaf->m_R - this_leaf->m_L)){
            int right_leaf_size = ((this_leaf->m_right_sibling)->m_R - (this_leaf->m_right_sibling)->m_L);
            Node* right_leaf = this_leaf->m_right_sibling;
            int last_filled_spot_in_arr = m_rightest_leaf->m_R;

            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[right_leaf->m_L + right_leaf->m_holes], NUMDIMS*(right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[right_leaf->m_L + right_leaf->m_holes], NUMDIMS*(right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + right_leaf->m_L + right_leaf->m_holes, (right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes)*sizeof(int));

            folan = (right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes);
            right_leaf->m_L = last_filled_spot_in_arr;
            right_leaf->m_R = last_filled_spot_in_arr + folan;

            memcpy(m_data_arr_mins[right_leaf->m_R], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[right_leaf->m_R], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + right_leaf->m_R, tba_ids + this_tba_start, k*sizeof(int));
            
            right_leaf->m_R += k;
            
            (right_leaf->m_left_sibling)->m_right_sibling = right_leaf->m_right_sibling;
            (right_leaf->m_right_sibling)->m_left_sibling = right_leaf->m_left_sibling;

            right_leaf->m_left_sibling = m_rightest_leaf;
            m_rightest_leaf->m_right_sibling = right_leaf;
            right_leaf->m_right_sibling = NULL;
            right_leaf->m_holes = 0;
            right_leaf->m_count = right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes;
            m_rightest_leaf = right_leaf;
            
            (this_leaf->m_right_sibling)->m_L -= right_leaf_size;
            (this_leaf->m_right_sibling)->m_holes += (right_leaf_size);
            return;
        }
        else{
            int this_L = this_leaf->m_L, this_R = this_leaf->m_R;
            ASSERT(this_L == 0);
            Node* right_leaf = this_leaf->m_right_sibling;
            (this_leaf->m_right_sibling)->m_left_sibling = NULL;
            m_leftest_leaf = right_leaf;

            int last_filled_spot_in_arr = m_rightest_leaf->m_R;

            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);

            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_L + this_leaf->m_holes], NUMDIMS*(folan)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_L + this_leaf->m_holes], NUMDIMS*(folan)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_L + this_leaf->m_holes, (folan)*sizeof(int));

            this_leaf->m_L = last_filled_spot_in_arr;
            this_leaf->m_R = last_filled_spot_in_arr + folan;

            this_leaf->m_right_sibling = NULL;
            this_leaf->m_left_sibling = m_rightest_leaf;
            m_rightest_leaf->m_right_sibling = this_leaf;
            this_leaf->m_holes = 0;
            m_rightest_leaf = this_leaf;
            
            memcpy(m_data_arr_mins[this_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_R - k, tba_ids + this_tba_start, k*sizeof(int));

            right_leaf->m_L =0;
            right_leaf->m_holes = this_R - k;
            ASSERT(right_leaf->m_id == m_leftest_leaf->m_id);
            return;
        }
    }

    else if(!is_leaf_in_OL(this_leaf, overlapping_leaves)){
        if(k > (this_leaf->m_R - this_leaf->m_L)){

            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));

            m_pending_insertions_count+=(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);

            memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, (this_leaf->m_R - this_leaf->m_L)*sizeof(int));
            (this_leaf->m_right_sibling)->m_L = this_leaf->m_L;

            Node *next_leaf = this_leaf->m_left_sibling; 
            int this_leaf_size = this_leaf->m_R - this_leaf->m_L;

            int branch_index = getNodeBranchIndex(this_leaf);
            Node* parent = this_leaf->m_parent;
            DisconnectBranch(parent, branch_index);
            deleteEmptyNodesUp(parent);
            
            shuffle_left(overlapping_leaves, next_leaf, this_tba_start + this_leaf_size, k - this_leaf_size);
        }
        else{
            int src = max(this_leaf->m_L + this_leaf->m_holes, this_leaf->m_R - this_leaf->m_holes);
            int how_many_to_move = min(this_leaf->m_holes,this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);

            memcpy(m_data_arr_mins[this_leaf->m_L], m_data_arr_mins[src], NUMDIMS*(how_many_to_move)*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], m_data_arr_maxes[src], NUMDIMS*(how_many_to_move)*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, m_data_arr_ids + src, how_many_to_move*sizeof(int));

            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_R - k], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_R - k], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_R - k, (k - this_leaf->m_holes)*sizeof(int));

            m_pending_insertions_count+=(k - this_leaf->m_holes);
            
            memcpy(m_data_arr_mins[this_leaf->m_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_R - k, tba_ids + this_tba_start, k*sizeof(int));
            
            this_leaf->m_R -= k;
            (this_leaf->m_right_sibling)->m_L -= k;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;

            (this_leaf->m_right_sibling)->m_holes = 0;
            return;
        }
    }
    else{
        if(k >= (this_leaf->m_R - this_leaf->m_L)){
            Node *next_leaf = this_leaf->m_left_sibling;
            (this_leaf->m_right_sibling)->m_L = this_leaf->m_L;
            (this_leaf->m_right_sibling)->m_holes = 0;
            
            int this_old_L = this_leaf->m_L, this_old_R = this_leaf->m_R;
            int last_filled_spot_in_arr = m_rightest_leaf->m_R;

            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));

            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
            this_leaf->m_L = last_filled_spot_in_arr;
            this_leaf->m_R = last_filled_spot_in_arr + folan;
            (this_leaf->m_left_sibling)->m_right_sibling = this_leaf->m_right_sibling;
            (this_leaf->m_right_sibling)->m_left_sibling = this_leaf->m_left_sibling;

            this_leaf->m_left_sibling = m_rightest_leaf;
            m_rightest_leaf->m_right_sibling = this_leaf;
            this_leaf->m_right_sibling = NULL;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            m_rightest_leaf = this_leaf;

            memcpy(m_data_arr_mins[this_old_L], tba_mins[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_old_L], tba_maxes[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_old_L, tba_ids + this_tba_start, (this_old_R - this_old_L)*sizeof(int));
            shuffle_left(overlapping_leaves, next_leaf, this_tba_start + this_old_R - this_old_L, k - (this_old_R - this_old_L));
        }
        else{
            int how_many_to_move = min(this_leaf->m_holes, this_leaf->m_R - this_leaf->m_holes-this_leaf->m_L);
            memcpy(m_data_arr_mins[this_leaf->m_L], m_data_arr_mins[this_leaf->m_R - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], m_data_arr_maxes[this_leaf->m_R - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, m_data_arr_ids + this_leaf->m_R - how_many_to_move, how_many_to_move*sizeof(int));

            float* temp_mins;
            temp_mins = (float*)malloc((k - this_leaf->m_holes) * NUMDIMS * sizeof(float));

            float* temp_maxes;
            temp_maxes = (float*)malloc((k - this_leaf->m_holes) * NUMDIMS * sizeof(float));

            int* temp_ids = (int*)malloc((k - this_leaf->m_holes) * sizeof(int));

            memcpy(temp_mins, m_data_arr_mins[this_leaf->m_R - k], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(temp_maxes, m_data_arr_maxes[this_leaf->m_R - k], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(temp_ids, m_data_arr_ids + this_leaf->m_R - k, (k - this_leaf->m_holes)*sizeof(int));

            
            memcpy(m_data_arr_mins[this_leaf->m_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_R - k, tba_ids + this_tba_start, k*sizeof(int));

            memcpy(tba_mins[this_tba_start], temp_mins, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(tba_maxes[this_tba_start], temp_maxes, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(tba_ids + this_tba_start, temp_ids, (k - this_leaf->m_holes)*sizeof(int));

            this_leaf->m_R -= k;
            (this_leaf->m_right_sibling)->m_L -= k;
            k -= this_leaf->m_holes;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            (this_leaf->m_right_sibling)->m_holes = 0;
            shuffle_left(overlapping_leaves, this_leaf->m_left_sibling, this_tba_start, k);

            free(temp_ids);
            free(temp_mins);
            free(temp_maxes);
        }
    }
}


RTREE_TEMPLATE
void RTREE_QUAL::ripple_v2(vector<overlapping_leaf_tbas> *overlapping_leaves){

    Node* this_leaf;
    int this_tba_count;
    int this_tba_start;
    bool dir_is_right;
    int k;
    int folan, filan;

    for(int index =0; index < overlapping_leaves->size(); index++){
        this_leaf = overlapping_leaves->at(index).this_leaf;
        this_tba_start = overlapping_leaves->at(index).this_tba_start;
        this_tba_count = overlapping_leaves->at(index).this_tba_count;
        k = this_tba_count;

        if(k == 0) continue;

        if(this_leaf->m_id == m_leftest_leaf->m_id){
            if(k <= this_leaf->m_holes){
                memcpy(m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes - k, tba_ids + this_tba_start, k*sizeof(int));
                memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
                this_leaf->m_holes -= k;
                this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
                continue;
            }
            // otherwose move it down
            // move leaf to end with tbi
            // continue
            Node* right_leaf = this_leaf->m_right_sibling;
            int old_R = this_leaf->m_R;
            int last_filled_spot = m_rightest_leaf->m_R;

            memcpy(m_data_arr_mins[last_filled_spot], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes] , NUMDIMS* (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes] , NUMDIMS* (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes , (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));

            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
            this_leaf->m_L = last_filled_spot;
            this_leaf->m_R = last_filled_spot + folan;

            memcpy(m_data_arr_ids + this_leaf->m_R, tba_ids + this_tba_start, k*sizeof(int));
            memcpy(m_data_arr_mins[this_leaf->m_R], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));

            this_leaf->m_R += k;

            this_leaf->m_left_sibling = m_rightest_leaf;
            this_leaf->m_right_sibling = NULL;
            m_rightest_leaf->m_right_sibling = this_leaf;
            m_rightest_leaf = this_leaf;
            this_leaf->m_holes = 0;

            right_leaf->m_left_sibling = NULL;
            m_leftest_leaf = right_leaf;
            right_leaf->m_L = 0;
            right_leaf->m_holes += old_R;

            continue;
            

        }
        else if(this_leaf->m_id == m_rightest_leaf->m_id){
            memcpy(m_data_arr_ids + this_leaf->m_R, tba_ids + this_tba_start, k*sizeof(int));
            memcpy(m_data_arr_mins[this_leaf->m_R], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));

            this_leaf->m_R += k;
            this_leaf->m_count += k;
            continue;
        }
        else{
            dir_is_right = true;

            if(this_leaf->m_id == m_leftest_leaf->m_id || is_leaf_in_OL(this_leaf->m_right_sibling, overlapping_leaves)){
                dir_is_right = false;
            }
            
            if(k < (this_leaf->m_R - this_leaf->m_L)){
                if(k <= this_leaf->m_holes){
                    memcpy(m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes - k, tba_ids + this_tba_start, k*sizeof(int));
                    memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
                    memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
                    this_leaf->m_holes -= k;
                    this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
                    continue;
                }
                memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], this_leaf->m_holes*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], this_leaf->m_holes*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, this_leaf->m_holes*sizeof(int));

                #ifdef timebreak
                    auto start_time = clock();
                #endif
                if(dir_is_right){
                    shuffle_right(overlapping_leaves, this_leaf->m_right_sibling, this_tba_start + this_leaf->m_holes, k - this_leaf->m_holes);
                    this_leaf->m_holes = 0;
                }
                else{
                    shuffle_left(overlapping_leaves, this_leaf->m_left_sibling, this_tba_start + this_leaf->m_holes, k - this_leaf->m_holes);
                }
                #ifdef timebreak
                    ripple_time += (clock() - start_time);
                #endif

                this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;


            }
            else{

                (this_leaf->m_right_sibling)->m_L = this_leaf->m_L;
                (this_leaf->m_right_sibling)->m_holes += (this_leaf->m_R - this_leaf->m_L);

                int last_filled_spot_in_arr = m_rightest_leaf->m_R;

                memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
                memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
                memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));

                folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
                this_leaf->m_L = last_filled_spot_in_arr;
                last_filled_spot_in_arr += folan;

                memcpy(m_data_arr_mins[last_filled_spot_in_arr], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_maxes[last_filled_spot_in_arr], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_ids + last_filled_spot_in_arr, tba_ids + this_tba_start, k*sizeof(int));

                this_leaf->m_R = last_filled_spot_in_arr + k;

                (this_leaf->m_left_sibling)->m_right_sibling = this_leaf->m_right_sibling;
                (this_leaf->m_right_sibling)->m_left_sibling = this_leaf->m_left_sibling;

                this_leaf->m_left_sibling = m_rightest_leaf;
                m_rightest_leaf->m_right_sibling = this_leaf;
                this_leaf->m_right_sibling = NULL;
                this_leaf->m_holes = 0;
                this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
                m_rightest_leaf = this_leaf;
            }
        }
    }
}


RTREE_TEMPLATE
bool RTREE_QUAL::CountDataPoints(Node *a_node, int &a_foundCount){
    if (a_node->m_parent == NULL) a_foundCount += m_pending_insertions_count;
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            CountDataPoints(a_node->m_branch[index].m_child, a_foundCount);
        }
    } else {
        a_foundCount += (a_node->m_R - a_node->m_L - a_node->m_holes);
    }
    return true; // Continue searching
}



RTREE_TEMPLATE
bool RTREE_QUAL::CountNodes(Node *a_node, int &a_foundCount){
    a_foundCount += 1;
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            CountNodes(a_node->m_branch[index].m_child, a_foundCount);
        }
    }
    return true; // Continue searching
}


RTREE_TEMPLATE
bool RTREE_QUAL::SumDataPointIDs(Node *a_node, int &a_foundCount){
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            SumDataPointIDs(a_node->m_branch[index].m_child, a_foundCount);
        }
    } else {
        for(int i= 0; i < a_node->m_count; i++){
            a_foundCount += a_node->m_branch[i].m_data;
        }
    }
    return true; // Continue searching
}



RTREE_TEMPLATE
int RTREE_QUAL::ChooseAxisLargestSideMid(Rect node_cover, Rect query_cover) {

    float this_ms;
    float current_max_m = std::abs(node_cover.m_max[0] - node_cover.m_min[0]);
    int current_max_index = 0;
    for(int i =1; i < NUMDIMS; i++){
        this_ms = std::abs(node_cover.m_max[i] - node_cover.m_min[i]);
        if(this_ms > current_max_m){
            current_max_m = this_ms;
            current_max_index = i;
        }
    }

    float mid_of_node_cover_axis = (float) ((node_cover.m_max[current_max_index] + node_cover.m_min[current_max_index]) / 2);
    if(std::abs(query_cover.m_min[current_max_index] - mid_of_node_cover_axis) <
    std::abs(query_cover.m_max[current_max_index] - mid_of_node_cover_axis)){
        return (2 * current_max_index);
    } else{
        return (2 * current_max_index + 1);
    }
}



RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMin_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v45(L, R, crack_value, axis, true, left_rect, right_rect);
}


RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMax_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v45(L, R, crack_value, axis, false, left_rect, right_rect);
}



RTREE_TEMPLATE
void RTREE_QUAL::swap(float &a, float &b) {
    float t = a;
    a = b;
    b = t;
}

RTREE_TEMPLATE
void RTREE_QUAL::swap_index(int i, int j) {
    std::swap(m_data_arr_ids[i], m_data_arr_ids[j]);
    std::swap(m_data_arr_mins[i], m_data_arr_mins[j]);
    std::swap(m_data_arr_maxes[i], m_data_arr_maxes[j]);
}




RTREE_TEMPLATE
int RTREE_QUAL::
MyPartition_v45(int low, int high, float pivot_value, int axis, bool min_or_max, Rect *left_rect, Rect *right_rect){
    int x1 = low; int x2 = high - 1;
    Rect x1_rect, x2_rect;

    if(min_or_max){
        while(x1 <= x2 && x2 > 0){
            if (m_data_arr_mins[x1][axis] <= pivot_value) {
                // update left cover

                if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];}

                x1++;
            }
            else{
                while(x2 > 0 && x2 >= x1 && m_data_arr_mins[x2][axis] >= pivot_value){
                    // update right cover

                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){right_rect->m_min[axis] = m_data_arr_mins[x2][axis];}
                    if(m_data_arr_maxes[x2][axis] > right_rect->m_max[axis]){right_rect->m_max[axis] = m_data_arr_maxes[x2][axis];}

                    x2--;
                }

                if(x1 < x2){
                    swap_index(x1, x2);
                    // update left cover

                    if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];}

                    // update right cover

                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){right_rect->m_min[axis] = m_data_arr_mins[x2][axis];}
                    if(m_data_arr_maxes[x2][axis] > right_rect->m_max[axis]){right_rect->m_max[axis] = m_data_arr_maxes[x2][axis];}

                    x2--; x1++;
                }
            }
        }
    }
    else{
        while(x1 <= x2 && x2 > 0){
            if (m_data_arr_maxes[x1][axis] < pivot_value) {
                // update left cover

                if(m_data_arr_mins[x1][axis] < left_rect->m_min[axis]){
                    left_rect->m_min[axis] = m_data_arr_mins[x1][axis];
                }
                if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){
                    left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];
                }

                x1++;

            }
            else{
                while(x2 > 0 && x2 >= x1 && m_data_arr_maxes[x2][axis] >= pivot_value){
                    // update right cover
                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){
                        right_rect->m_min[axis] = m_data_arr_mins[x2][axis];
                    }

                    x2--;
                }
                if(x1 < x2){
                    swap_index(x1, x2);
                    // update left cover

                    if(m_data_arr_mins[x1][axis] < left_rect->m_min[axis]){
                        left_rect->m_min[axis] = m_data_arr_mins[x1][axis];
                    }
                    if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){
                        left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];
                    }

                    // update right cover
                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){
                        right_rect->m_min[axis] = m_data_arr_mins[x2][axis];
                    }
                    x2--; x1++;
                }
            }
        }
    }

    return x1;
}


RTREE_TEMPLATE
void RTREE_QUAL::CopyBranch(Branch &current, const Branch &other){
    current.m_rect = other.m_rect;
    current.m_child = other.m_child;
    current.m_data = other.m_data;
}


RTREE_TEMPLATE
int RTREE_QUAL::getNodeBranchIndex(Node *a_node) {
    for(int i = 0; i < (a_node->m_parent)->m_count; i++) {
        if (a_node->m_id == ((a_node->m_parent)->m_branch[i].m_child)->m_id) {
            return i;
        }
    }
    return -1;
}


RTREE_TEMPLATE
void RTREE_QUAL::Insert_anylevel(const Branch &a_branch, Node *start, int a_level){
    ASSERT(a_branch.m_child->m_count > 0);
    if(start->m_level > a_level){
        int index = PickBranch(&a_branch.m_rect, start);

        Insert_anylevel(a_branch, start->m_branch[index].m_child, a_level);
    }
    else{
        Node * current_node = start;
        Branch current_branch;
        CopyBranch(current_branch, a_branch);
        bool is_split = false; // so the loop starts
        bool just_go_till_top = false; // to indicate that no more adds are required, but we should recurse up to update the rects
        while(true){
            ASSERT(current_branch.m_child->m_count > 0);
            if(!just_go_till_top) {
                Node *newNode = AllocNode();
                Rect rect1; Rect rect2;
                is_split = AddBranch(&current_branch, current_node, &newNode, rect1, rect2);

                if (is_split) {
                    // was split
                    if (current_node->m_parent != NULL) {
                        ASSERT(newNode->m_count > 0);
                        
                        int index = getNodeBranchIndex(current_node);
                        (current_node->m_parent)->m_branch[index].m_rect = rect1;
                        current_branch.m_child = newNode;
                        
                        current_branch.m_rect = rect2;
                        current_node = current_node->m_parent;
                    } else {
                        // now we have to make a new root

                        Node *newRoot = AllocNode();
                        newRoot->m_level = current_node->m_level + 1;
                        newRoot->m_parent = NULL;

                        current_branch.m_child = newNode;
                        current_branch.m_rect = rect2;

                        AddBranch(&current_branch, newRoot, NULL);

                        current_branch.m_child = current_node;
                        current_branch.m_rect = rect1;

                        AddBranch(&current_branch, newRoot, NULL);

                        m_root = newRoot;
                        break;
                    }
                } else {
                    // was not split
                    if (current_node->m_parent != NULL) {
                        int index = getNodeBranchIndex(current_node);
                        (current_node->m_parent)->m_branch[index].m_rect = CombineRect(&(current_branch.m_rect),
                                                                                       &((current_node->m_parent)->m_branch[index].m_rect));
                        current_node = current_node->m_parent;
                        just_go_till_top = true;
                    } else {
                        break;
                    }
                }
            } else{
                if (current_node->m_parent != NULL) {
                    int index = getNodeBranchIndex(current_node);
                    (current_node->m_parent)->m_branch[index].m_rect = CombineRect(&(current_branch.m_rect),
                                                                                   &((current_node->m_parent)->m_branch[index].m_rect));
                    current_node = current_node->m_parent;
                } else{
                    break;
                }
            }

        }

    }
}



RTREE_TEMPLATE
void RTREE_QUAL::PrintLeafSizesRec(Node* a_node, ofstream &myfile){
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            PrintLeafSizesRec(a_node->m_branch[index].m_child, myfile);
        }
    } else{
        myfile << a_node->m_id << " " << a_node->m_count << "\n";
    }

}


RTREE_TEMPLATE
void RTREE_QUAL::printLeafLinkedList(){
    cout << "LEAF LINKED LIST ---------------" << endl;
    Node *current_node = m_leftest_leaf;
    while(current_node != NULL){
        cout << "node " << current_node->m_id << " L: " << current_node->m_L << " R: " << current_node->m_R << " holes: " << current_node->m_holes << " count: " << current_node->m_count << endl;
        current_node = current_node->m_right_sibling;
    }
    cout << "END of LEAF LINKED LIST ---------------" << endl;
}


RTREE_TEMPLATE
void RTREE_QUAL::printIDs(){
    cout << "pending insertion ids: ... " << endl;
    for(int i = 0; i < m_pending_insertions_count; i++){
        cout << m_pending_insertions_ids[i] << " " ;
    }
    cout << endl << "----------" << endl;
    cout << "data arr ids: ... " << endl;
    for(int i = 0; i < 2 * DATA_COUNT; i++){
        cout << m_data_arr_ids[i] << " ";
    }
    cout << endl << "----------" << endl;

}


RTREE_TEMPLATE
void RTREE_QUAL::printLeafBoundsOfData(int data_index){
    printLeafBoundsOfDataRec(m_root, data_index);
}

RTREE_TEMPLATE
bool RTREE_QUAL::printLeafBoundsOfDataRec(Node *a_node, int data_index){
    if(a_node->IsInternalNode()){
        for(int i = 0; i < a_node->m_count; i++){
            bool print_or_not = printLeafBoundsOfDataRec(a_node->m_branch[i].m_child, data_index);
            if(print_or_not){
                cout << "level: " << a_node->m_level << " : " ;
                cout << "child's level: " << a_node->m_branch[i].m_child->m_level << " : ";
                for(int j = 0; j < NUMDIMS; j++){
                    cout << a_node->m_branch[i].m_rect.m_min[j] << " " << a_node->m_branch[i].m_rect.m_max[j] << " , ";
                }
                cout << endl;
                return true;
            }
        }
    }
    else{
        for(int j = a_node->m_L + a_node->m_holes; j < a_node->m_R; j++){
            if(m_data_arr_ids[j] == data_index){
                cout << "data was: " << endl;
                for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[j][jj] << " " << m_data_arr_maxes[j][jj] << " ,";
                cout << endl; 
                return true;
            }
        }
        return false;
    }
    return false;
}


RTREE_TEMPLATE
void RTREE_QUAL::findAllLeavesWithDataInThem(int data_index){
    Node* this_leaf = m_leftest_leaf;
    while(this_leaf != NULL){
        for(int i = this_leaf->m_L + this_leaf->m_holes; i < this_leaf->m_R; i++){
            if(m_data_arr_ids[i] == data_index){
                cout << "data in leaf " << this_leaf->m_id << " is item " << i << " and L, holes, R: " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;
                cout << "data bounds: " << m_data_arr_mins[i][0] << " " << m_data_arr_maxes[i][0] << ", " << m_data_arr_mins[i][1] << " " << m_data_arr_maxes[i][1] << endl;
                if(this_leaf->m_parent != NULL){
                    int branch_index = getNodeBranchIndex(this_leaf);
                    Rect cover = (this_leaf->m_parent)->m_branch[branch_index].m_rect;
                    cout << "leaf bound: " << cover.m_min[0] << " " << cover.m_min[1] << " " << cover.m_max[0] << " " << cover.m_max[1] << endl;
                }
            }
        }
        // move on to the next
        this_leaf = this_leaf->m_right_sibling;
    }

    // also look in pending
    for(int i = 0; i < m_pending_insertions_count; i++){
        if(m_pending_insertions_ids[i] == data_index){
            cout << "data is object " << i << " in pending" << endl;
            cout << "data bounds: " << m_pending_insertions_mins[i][0] << " " << m_pending_insertions_maxes[i][0] << ", " << m_pending_insertions_mins[i][1] << " " << m_pending_insertions_maxes[i][1] << endl;

        }
    }
    cout << "FINISHED\n";
}

RTREE_TEMPLATE
void RTREE_QUAL::lookForEmptyNodesRec(Node *a_node){

    if(a_node->IsInternalNode()){
        if(a_node->m_count > 0){
            for(int i = 0; i < a_node->m_count; i++){
                    lookForEmptyNodesRec((a_node->m_branch[i]).m_child);
            }
        }
        else{
            cout << "node " << a_node->m_id << " was empty. very very bad. not recommended. BOOHOO" << endl;
            cout << "info: count = " << a_node->m_count << " level: " << a_node->m_level << " MAX-lvl " << m_root->m_level << endl;
            // return;
            exit(3);
        }
    }
    return;
}

RTREE_TEMPLATE
void RTREE_QUAL::lookForEmptyNodes(){
    lookForEmptyNodesRec(m_root);
}

RTREE_TEMPLATE
void RTREE_QUAL::findNodeCount(int node_id, Node* this_node){
    if(this_node->m_id == node_id){
        cout << "found node, has " << this_node->m_count << " count, lvl: " << this_node->m_level << endl;
        return;
    }
    if(this_node->IsInternalNode()){
        if(this_node->m_count > 0){
            for(int i = 0; i < this_node->m_count; i++){
                    findNodeCount(node_id, (this_node->m_branch[i]).m_child);
            }
        }
    }
}


RTREE_TEMPLATE
std::pair<int, int> RTREE_QUAL::CountLeaves(){
    int reg_counter = 0; int irreg_counter = 0;
    Node* this_leaf = m_leftest_leaf;
    while(this_leaf != NULL){
        if(this_leaf->isRegular())  reg_counter++;
        else irreg_counter++;
        // go to next leaf
        this_leaf = this_leaf->m_right_sibling;
    }
    return std::make_pair(reg_counter, irreg_counter);
}

#undef RTREE_TEMPLATE
#undef RTREE_QUAL

#endif //GLIDE_CQRIPPLE_H

