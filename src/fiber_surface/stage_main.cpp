#include "pipeline.h"

// Each target selects one independently runnable stage. The actual algorithms
// are compiled once from their five named .cpp files into fiber_surface_core.
int main(int argc,char** argv) {
    return fiber_surface::stage_main(static_cast<fiber_surface::Stage>(FIBER_STAGE),argc,argv);
}
