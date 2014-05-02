/**
   rtree lib usage example app.
*/

#include <stdio.h>
#include "rtree.h"

RTREEMBR rects[] = {
    { {0, 0, 0, 2, 2, 0} },  /* xmin, ymin, zmin, xmax, ymax, zmax (for 3 dimensional RTree) */
    { {5, 5, 0, 7, 7, 0} },
    { {8, 5, 0, 9, 6, 0} },
    { {7, 1, 0, 9, 2, 0} }
};


int nrects = sizeof(rects) / sizeof(rects[0]);
RTREEMBR search_rect = {
    {6, 4, 0, 10, 6, 0}   /* search will find above rects that this one overlaps */
};

int MySearchCallback(int id, void* arg) 
{
    /* Note: -1 to make up for the +1 when data was inserted */
    fprintf (stdout, "Hit data mbr %d ", id-1);
    return 1; /* keep going */
}

int main()
{
    RTREENODE* root = RTreeCreate();
    
    int i, nhits;
    
    fprintf (stdout, "nrects = %d ", nrects);
    
    /* Insert all the testing data rects */
    for(i=0; i<nrects; i++){
        RTreeInsertRect(&rects[i],  /* the mbr being inserted */
                        i+10,        /* i+1 is mbr ID. ID MUST NEVER BE ZERO */
                        &root,        /* the address of rtree's root since root can change undernieth*/
                        0            /* always zero which means to add from the root */
            );
    }

    nhits = RTreeSearch(root, &search_rect, MySearchCallback, 0);
    
    fprintf (stdout, "Search resulted in %d hits ", nhits);

    RTreeDestroy (root);

    return 0;
}
