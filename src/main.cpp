
#include "b_ism.h"
#include "compare_results.h"
#include "exec_timer.h"

int main_wo_r() {

// Pre simulated Data given the above defaults
    Eigen::MatrixXd data(3 * 10, 6);

    data << 1, 1, -21.60672924, -2.30936055, -0.55970183, 1.34953562,
            1, 2, -13.53344759, -1.46139144, 0.21277444, -0.47043324,
            1, 3, 2.89082517, -0.06862945, -0.12904591, 1.80418798,
            1, 4, 2.62529980, -0.23848014, 0.31444212, 0.18225964,
            1, 5, -0.23413276, 0.54917267, 0.18824194, -2.76781820,
            1, 6, 7.73232752, 1.04541866, 0.01999001, -1.02430597,
            1, 7, -3.82418503, -0.37946788, 0.32571192, -1.17450487,
            1, 8, 4.97869714, -0.22696702, -0.28827137, 1.35956658,
            1, 9, -5.35322334, 0.09199853, -0.72593770, -1.00533029,
            1, 10, -0.68939691, -0.85342818, -0.29560190, 1.31471504,
            2, 1, -14.61152255, -0.50768911, -0.93493273, -0.60211238,
            2, 2, 8.11772892, 1.12089459, -0.43084369, 0.55911127,
            2, 3, -3.01870024, 0.58630229, -0.68998967, -0.40238799,
            2, 4, -13.25864790, -1.27681905, 0.06944209, -0.88038903,
            2, 5, 3.44635506, -0.70259846, 0.75558385, 0.55643298,
            2, 6, 25.48550732, 1.69326930, 1.02409623, -0.06464430,
            2, 7, -14.27278142, -0.27937351, -0.68313724, -1.14558556,
            2, 8, 16.07181146, 1.23860288, 0.89303727, -1.17465964,
            2, 9, -19.91590743, -0.27851469, -1.51865248, -0.30728423,
            2, 10, 27.03493229, 0.50750409, 2.13166805, 0.07164625,
            3, 1, -14.29968334, -0.52964157, -0.84032731, -0.24943370,
            3, 2, 29.52125002, 1.24425732, 1.89566249, -0.07214188,
            3, 3, 5.87508461, 0.74228208, 0.57532928, -2.24319079,
            3, 4, 2.27248084, -0.46578083, 0.89250538, -0.64141230,
            3, 5, 14.90452209, 0.75886611, 0.61983879, 0.82792110,
            3, 6, 7.14433418, 0.45025558, 0.09238917, 0.81542852,
            3, 7, 0.90919475, 0.61893994, -0.07080554, -1.25784309,
            3, 8, 5.77267567, 0.43634622, -0.35005450, -0.11806201,
            3, 9, 1.53707900, -0.06634003, -0.35659074, -0.11128336,
            3, 10, 2.72656784, -0.24980976, -0.11055986, 0.75599906;


    bool include_constant = true;
    bool tfe = false;
    bool ife = false;
    bool iis = true;
    bool sis = true;
    int i_index = 0; // -1 from R code because of 0 Indexing
    int t_index = 1; // -1 from R code because of 0 Indexing
    int y_index = 2; // -1 from R code because of 0 Indexing
    long Ndraw = 1000L;
    long Nburn = 500L;
    double lambda_b = 1000;
    double c0 = 0.0001;
    double C0 = 0.0001;

    bool geweke = false;
    std::string b_prior = "g";

    BismResults result1;
    BismResults result2;

    FunctionTimer function_timer = FunctionTimer();

//    function_timer.register_and_run("No optimization", b_ism,
//                                    data,
//                                    include_constant,
//                                    tfe,
//                                    ife,
//                                    iis,
//                                    sis,
//                                    y_index,
//                                    i_index,
//                                    t_index,
//                                    Ndraw,
//                                    Nburn,
//                                    b_prior,
//                                    lambda_b,
//                                    c0,
//                                    C0,
//                                    geweke,
//                                    result1,
//                                    ModelSelectionVersion::NO_OPTIMIZATION
//    );

    function_timer.register_and_run("Parallel Z", b_ism,
                                    data,
                                    include_constant,
                                    tfe,
                                    ife,
                                    iis,
                                    sis,
                                    y_index,
                                    i_index,
                                    t_index,
                                    Ndraw,
                                    Nburn,
                                    b_prior,
                                    lambda_b,
                                    c0,
                                    C0,
                                    geweke,
                                    result1,
                                    ModelSelectionVersion::SPLIT_MATRIX_PARALLEL
    );

    function_timer.register_and_run("Split z", b_ism,
                                    data,
                                    include_constant,
                                    tfe,
                                    ife,
                                    iis,
                                    sis,
                                    y_index,
                                    i_index,
                                    t_index,
                                    Ndraw,
                                    Nburn,
                                    b_prior,
                                    lambda_b,
                                    c0,
                                    C0,
                                    geweke,
                                    result2,
                                    ModelSelectionVersion::SPLIT_MATRIX
    );

    function_timer.print_function_summary();

    compare_bism_results(result1, result2, "No Optimization", "Split Z");

    return 0;
}

int main([[maybe_unused]] int argc, char *argv[]) {

    main_wo_r();
    return 0;
}

