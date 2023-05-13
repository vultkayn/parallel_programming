#ifndef _coordinate_h
#define _coordinate_h

#include <stdbool.h>

#ifndef __cplusplus
struct cord {
    float x0 ;
    float x1 ;
    float y0 ;
    float y1 ;
} ;

struct part_cord {
    float x ;
    float y ;
    float vx ;
    float vy ;
    char collide ;
} ;
#else
struct cord {
    float x0 = 0.0;
    float x1 = 0.0;
    float y0 = 0.0;
    float y1 = 0.0;
} ;

struct part_cord {
    float x = 0.0;
    float y = 0.0;
    float vx = 0.0;
    float vy = 0.0;
    char collide = 0;
} ;
#endif

typedef struct cord cord_t ;
typedef struct part_cord pcord_t ;

#endif
