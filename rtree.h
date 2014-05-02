/****************************************************************************
* RTree.H
*
* MODULE:       R-Tree library 
*              
* AUTHOR(S):    Antonin Guttman - original code
*               Daniel Green (green@superliminal.com) - major clean-up
*                               and implementation of bounding spheres
*               
* PURPOSE:      Multi Dimensional Index
*
* COPYRIGHT:    (C) 2001 by the GRASS Development Team
*
*               This program is free software under the GNU General Public
*               License (>=v2). Read the file COPYING that comes with GRASS
*               for details.
*
* LAST MODIFY:         ZhangLiang (cheungmine@gmail.com) - 2007-11
*****************************************************************************/
#ifndef  RTREE_H_INCLUDED
#define  RTREE_H_INCLUDED

/* PAGE_SIZE is normally the natural page size of the machine */
#define  PAGE_SIZE    512
#define  DIMS_NUMB    3       /* number of dimensions */
#define  SIDES_NUMB   2*DIMS_NUMB

/* typedef float REALTYPE; */
typedef double REALTYPE;


#ifndef  TRUE
#define  TRUE        1
#define  FALSE        0
#endif


typedef struct _RTREEMBR
{
    REALTYPE bound[SIDES_NUMB]; /* xmin,ymin,...,xmax,ymax,... */
}RTREEMBR;

typedef struct _RTREEBRANCH
{        
    RTREEMBR    mbr;
    struct _RTREENODE *child;    /* mbr id */
}RTREEBRANCH;

/* max branching factor of a node */
#define MAXCARD (int)((PAGE_SIZE-(2*sizeof(int))) / sizeof(RTREEBRANCH))

typedef struct _RTREENODE
{
    int    count;
    int    level; /* 0 is leaf, others positive */
    RTREEBRANCH  branch[MAXCARD];
}RTREENODE;

typedef struct _RTREELISTNODE
{
     struct _RTREELISTNODE    *next;
     RTREENODE        *node;
}RTREELISTNODE;

/*
* If passed to a tree search, this callback function will be called
* with the ID of each data mbr that overlaps the search mbr
* plus whatever user specific pointer was passed to the search.
* It can terminate the search early by returning 0 in which case
* the search will return the number of hits found up to that point.
*/
typedef int (*pfnSearchHitCallback)(int id, void* pfnParam);


int RTreeSetNodeMax(int new_max);

int RTreeSetLeafMax(int new_max);

int RTreeGetNodeMax(void);

int RTreeGetLeafMax(void);

/**
 * Initialize a rectangle to have all 0 coordinates.
 */
void RTreeInitRect( RTREEMBR *rc);

/**
 * Return a mbr whose first low side is higher than its opposite side -
 * interpreted as an undefined mbr.
 */
RTREEMBR RTreeNullRect(void);


/**
 * Print out the data for a rectangle.
 */
void RTreePrintRect( RTREEMBR *rc, int depth );

/**
 * Calculate the 2-dimensional area of a rectangle
 */
REALTYPE RTreeRectArea( RTREEMBR *rc );

/**
 * Calculate the n-dimensional volume of a rectangle
 */
REALTYPE RTreeRectVolume( RTREEMBR *rc );


/**
 * Calculate the n-dimensional volume of the bounding sphere of a rectangle
 * The exact volume of the bounding sphere for the given RTREEMBR.
 */
REALTYPE RTreeRectSphericalVolume( RTREEMBR *rc );


/**
 * Calculate the n-dimensional surface area of a rectangle
 */
REALTYPE RTreeRectSurfaceArea( RTREEMBR *rc );


/**
 * Combine two rectangles, make one that includes both.
 */
RTREEMBR RTreeCombineRect( RTREEMBR *rc1, RTREEMBR *rc2 );


/**
 * Decide whether two rectangles overlap.
 */
int RTreeOverlap( RTREEMBR *rc1, RTREEMBR *rc2);


/**
 * Decide whether rectangle r is contained in rectangle s.
 */
int RTreeContained( RTREEMBR *r, RTREEMBR *s);

/**
 * Split a node.
 * Divides the nodes branches and the extra one between two nodes.
 * Old node is one of the new ones, and one really new one is created.
 * Tries more than one method for choosing a partition, uses best result.
 */
void RTreeSplitNode( RTREENODE *node, RTREEBRANCH *br, RTREENODE **new_node);

/**
 * Initialize a RTREENODE structure. 
 */
void RTreeInitNode( RTREENODE *node );

/** 
 * Make a new node and initialize to have all branch cells empty. 
 */
RTREENODE *RTreeNewNode(void);

void RTreeFreeNode( RTREENODE *node );


/**
 * Print out the data in a node. 
 */
void RTreePrintNode( RTREENODE *node, int depth );


/**
 * Find the smallest rectangle that includes all rectangles in branches of a node.
 */
RTREEMBR RTreeNodeCover( RTREENODE *node );


/**
 * Pick a branch.  Pick the one that will need the smallest increase
 * in area to accomodate the new rectangle.  This will result in the
 * least total area for the covering rectangles in the current node.
 * In case of a tie, pick the one which was smaller before, to get
 * the best resolution when searching.
 */
int RTreePickBranch( RTREEMBR *rc, RTREENODE *node);


/**
 * Add a branch to a node.  Split the node if necessary.
 * Returns 0 if node not split.  Old node updated.
 * Returns 1 if node split, sets *new_node to address of new node.
 * Old node updated, becomes one of two.
 */
int RTreeAddBranch( RTREEBRANCH *br, RTREENODE *node, RTREENODE **new_node);


/**
 * Disconnect a dependent node. 
 */
void RTreeDisconnectBranch( RTREENODE *node, int i );


/**
 * Destroy (free) node recursively. 
 */
void RTreeDestroyNode ( RTREENODE *node );


/**
 * Create a new rtree index, empty. Consists of a single node. 
 */
RTREENODE * RTreeCreate(void);


/**
 * Destroy a rtree root must be a root of rtree. Free all memory.
 */
void RTreeDestroy(RTREENODE *root);


/**
 * Search in an index tree or subtree for all data rectangles that overlap the argument rectangle.
 * Return the number of qualifying data rects.
 */
int RTreeSearch( RTREENODE *node, RTREEMBR *rc, pfnSearchHitCallback pfnSHCB, void* pfnParam);

/** 
 * Insert a data rectangle into an index structure.
 * RTreeInsertRect provides for splitting the root;
 * returns 1 if root was split, 0 if it was not.
 * The level argument specifies the number of steps up from the leaf
 * level to insert; e.g. a data rectangle goes in at level = 0.
 * _RTreeInsertRect does the recursion.
 */
int RTreeInsertRect( RTREEMBR *rc, int tid, RTREENODE **root, int level);

/**
 * Delete a data rectangle from an index structure.
 * Pass in a pointer to a RTREEMBR, the tid of the record, ptr to ptr to root node.
 * Returns 1 if record not found, 0 if success.
 * RTreeDeleteRect provides for eliminating the root.
 */
int RTreeDeleteRect( RTREEMBR *rc, int tid, RTREENODE **root);

#endif /* RTREE_H_INCLUDED */

