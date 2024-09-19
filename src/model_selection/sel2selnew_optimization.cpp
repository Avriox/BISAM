//
// Created by jakob on 9/14/24.
//

#include "sel2selnew_optimization.h"

void
sel2selnew(int newgroup, int *sel, int *nsel, int *selnew, int *nselnew, bool copylast, int *ngroups, int *nvaringroup,
           int *firstingroup) {
//Copy sel into selnew. Elements are always kept ordered so that sel[0] < sel[1] < sel[2] ...
// - If newgroup in sel, don't copy it in selnew and set nselnew=nsel-1.
// - If newgroup not in sel, add it to selnew and set nselnew=nsel+1.
// - If copylast==true, copy last element sel[nsel] into selnew[nselnew]
    bool found;
    int i, ii, iii;
    for (i = 0, found = false; (i < *nsel) && (sel[i] <= firstingroup[newgroup]) && (!found); i++) {
        selnew[i] = sel[i];
        found = (sel[i] == firstingroup[newgroup]);
    }
    if (!found) { //add new group
        for (ii = 0; ii < nvaringroup[newgroup]; ii++) { selnew[i + ii] = firstingroup[newgroup] + ii; }
        for (iii = 0; (i + iii) < *nsel; iii++) { selnew[i + ii + iii] = sel[i + iii]; }
        (*nselnew) = (*nsel) + nvaringroup[newgroup];
    } else {  //remove new elem
        for (ii = i - 1 + nvaringroup[newgroup]; ii < *nsel; ii++) { selnew[ii - nvaringroup[newgroup]] = sel[ii]; }
        (*nselnew) = (*nsel) - nvaringroup[newgroup];
    }
    if (copylast) selnew[*nselnew] = sel[*nsel];
}

void sel2selnew_small_matrix(int newgroup, int *sel, int *nsel, int *selnew, int *nselnew, bool copylast, int *ngroups,
                             int *nvaringroup,
                             int *firstingroup) {
//Copy sel into selnew. Elements are always kept ordered so that sel[0] < sel[1] < sel[2] ...
// - If newgroup in sel, don't copy it in selnew and set nselnew=nsel-1.
// - If newgroup not in sel, add it to selnew and set nselnew=nsel+1.
// - If copylast==true, copy last element sel[nsel] into selnew[nselnew]
    bool found;
    int i, ii, iii;
    for (i = 0, found = false; (i < *nsel) && (sel[i] <= firstingroup[newgroup]) && (!found); i++) {
        selnew[i] = sel[i];
        found = (sel[i] == firstingroup[newgroup]);
    }
    if (!found) { //add new group
        for (ii = 0; ii < nvaringroup[newgroup]; ii++) {
            selnew[i + ii] = firstingroup[newgroup] + ii;
        }
        for (iii = 0; (i + iii) < *nsel; iii++) {
            selnew[i + ii + iii] = sel[i + iii];
        }
        (*nselnew) = (*nsel) + nvaringroup[newgroup];
    } else {  //remove new elem
        for (ii = i - 1 + nvaringroup[newgroup]; ii < *nsel; ii++) {
            selnew[ii - nvaringroup[newgroup]] = sel[ii];
        }
        (*nselnew) = (*nsel) - nvaringroup[newgroup];
    }
    if (copylast) selnew[*nselnew] = sel[*nsel];
}