//
// Created by jakob on 9/14/24.
//

#ifndef CPP_SEL2SELNEW_OPTIMIZATION_H
#define CPP_SEL2SELNEW_OPTIMIZATION_H

// Baseline implemenation from mombf
void sel2selnew(int newgroup, int *sel, int *nsel, int *selnew, int *nselnew, bool copylast, int *ngroups,
                int *nvaringroup, int *firstingroup);

// Optimized to break up the large matrix into smaller sub matrices and invert them
// individually
void sel2selnew_small_matrix(int newgroup, int *sel, int *nsel, int *selnew, int *nselnew, bool copylast, int *ngroups,
                             int *nvaringroup, int *firstingroup);

#endif //CPP_SEL2SELNEW_OPTIMIZATION_H
