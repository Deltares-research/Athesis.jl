// This file is a part of Julia. License is MIT: https://julialang.org/license

#include <julia.h>
#include <stdio.h>
#include <math.h>


#if defined(_WIN32)
#if !defined(MKERNEL_API)
#define ATHESIS_API __declspec(dllexport)
#endif
#else
#define ATHESIS_API __attribute__((visibility("default")))
#endif

jl_value_t *checked_eval_string(const char* code)
{
    jl_value_t *result = jl_eval_string(code);
    if (jl_exception_occurred()) {
        // none of these allocate, so a gc-root (JL_GC_PUSH) is not necessary
        jl_call2(jl_get_function(jl_base_module, "showerror"),
                 jl_stderr_obj(),
                 jl_exception_occurred());
        jl_printf(jl_stderr_stream(), "\n");
        jl_atexit_hook(1);
        exit(1);
    }
    assert(result && "Missing return value but no exception occurred!");
    return result;
}

ATHESIS_API int athesis_bmi_initialize()
{
    jl_init();
    checked_eval_string("import Pkg");
    checked_eval_string("Pkg.activate(\"..\")");
    checked_eval_string("Pkg.instantiate()");
    checked_eval_string("using Athesis");
    checked_eval_string("using BasicModelInterface");
    checked_eval_string("const BMI = BasicModelInterface");
    checked_eval_string("simulation = BMI.initialize(Simulation, nothing)");
    
}

ATHESIS_API int athesis_bmi_update()
{
    checked_eval_string("BMI.update(simulation)");
}

// int main()
// {
//     jl_init();

//     {
//         // Test Athesis
//         checked_eval_string("import Pkg");
//         checked_eval_string("Pkg.activate(\"..\")");
//         checked_eval_string("Pkg.instantiate()");
//         checked_eval_string("using Athesis");
//         checked_eval_string("using BasicModelInterface");
//         checked_eval_string("const BMI = BasicModelInterface");
//         checked_eval_string("simulation = BMI.initialize(Simulation, nothing)");
//         checked_eval_string("BMI.update(simulation)");
//     }

//     int ret = 0;
//     jl_atexit_hook(ret);
//     return ret;
// }
