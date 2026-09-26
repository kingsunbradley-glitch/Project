// macOS-specific entry point.  The shared implementation remains in the
// original macro so future physics/layout edits do not have to be duplicated.
#include "Draw_v4_all_combined.C"

void Draw_v4_all_combined_mac()
{
    Draw_v4_all_combined_impl(true);
}
