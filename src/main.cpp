
#include "./b_ism.h"
#include "./compare_results.h"
#include "./exec_timer.h"

int main_wo_r() {

// Pre simulated Data given the above defaults
    Eigen::MatrixXd data(2 * 10, 6);

    data << 1, 1, -16.29140398, -2.30936055, -0.55970183, 1.34953562,
            1, 2, -8.67599581, -1.46139144, 0.21277444, -0.47043324,
            1, 3, -2.78436407, -0.06862945, -0.12904591, 1.80418798,
            1, 4, -0.98097547, -0.23848014, 0.31444212, 0.18225964,
            1, 5, 5.89555601, 0.54917267, 0.18824194, -2.76781820,
            1, 6, 7.58003256, 1.04541866, 0.01999001, -1.02430597,
            1, 7, -1.22913630, -0.37946788, 0.32571192, -1.17450487,
            1, 8, 1.64000650, -0.22696702, -0.28827137, 1.35956658,
            1, 9, 4.05506831, 0.09199853, -0.72593770, -1.00533029,
            1, 10, -1.53000515, -0.85342818, -0.29560190, 1.31471504,
            2, 1, -2.84301984, -0.50768911, -0.93493273, -0.60211238,
            2, 2, 7.23214651, 1.12089459, -0.43084369, 0.55911127,
            2, 3, 3.24958704, 0.58630229, -0.68998967, -0.40238799,
            2, 4, -6.65240808, -1.27681905, 0.06944209, -0.88038903,
            2, 5, -3.58953437, -0.70259846, 0.75558385, 0.55643298,
            2, 6, 11.50349176, 1.69326930, 1.02409623, -0.06464430,
            2, 7, 0.41077442, -0.27937351, -0.68313724, -1.14558556,
            2, 8, 9.39463975, 1.23860288, 0.89303727, -1.17465964,
            2, 9, -0.86074003, -0.27851469, -1.51865248, -0.30728423,
            2, 10, 1.36886443, 0.50750409, 2.13166805, 0.07164625;


    bool include_constant = false;
    bool tfe = true;
    bool ife = true;
    bool iis = true;
    bool sis = true;
    int i_index = 0; // -1 from R code because of 0 Indexing
    int t_index = 1; // -1 from R code because of 0 Indexing
    int y_index = 2; // -1 from R code because of 0 Indexing
    long Ndraw = 5000L;
    long Nburn = 1000L;
    double lambda_b = 1000;
    double c0 = 0.0001;
    double C0 = 0.0001;

    double tau = 1;
    double va = 1;
    double vb = 1;

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
                                    tau,
                                    va,
                                    vb,
                                    ModelSelectionVersion::SPLIT_MATRIX
    );

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
                                    tau,
                                    va,
                                    vb,
                                    ModelSelectionVersion::SPLIT_MATRIX_PARALLEL
    );


    function_timer.print_function_summary();

    compare_bism_results(result1, result2, "No Optimization", "Split Z");

    return 0;
}

int main([[maybe_unused]] int argc, char *argv[]) {

    main_wo_r();
    return 0;
}

