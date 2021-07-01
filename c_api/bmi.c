// This file is a part of Julia. License is MIT: https://julialang.org/license

#include <julia.h>
#include <stdio.h>


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

void check_exception()
{
    if (jl_exception_occurred()) {
        // none of these allocate, so a gc-root (JL_GC_PUSH) is not necessary
        jl_call2(jl_get_function(jl_base_module, "showerror"),
                 jl_stderr_obj(),
                 jl_exception_occurred());
        jl_printf(jl_stderr_stream(), "\n");
        jl_atexit_hook(1);
        exit(1);
    }
}

ATHESIS_API int initialize(char* config_file)
{   
    // Initialize Julia
    jl_init();

    // Activate Athesis
    checked_eval_string("import Pkg");
    // The path to Athesis is hardcoded here
    // If the C API becomes more than a proof of concept at some point that should be changed
    checked_eval_string("Pkg.activate(\"C:/checkouts/Athesis/Athesis.jl\")");
    checked_eval_string("Pkg.instantiate()");
    checked_eval_string("using Athesis");
    checked_eval_string("using BasicModelInterface");
    checked_eval_string("const BMI = BasicModelInterface");

    // Initialize via the BMI
    char* command = (char*)malloc(256 * sizeof(char));
    sprintf(command, "simulation = BMI.initialize(Simulation, \"%s\")", config_file);
    checked_eval_string(command);
    free(command);
    return 0;
}

// Note that `update`should not contain any parameters according to BMI
// Since the goal is to couple to dimr, we follow this standard though
ATHESIS_API int update(double dt)
{
    // Rebinding method to global namespace, until I found out how to do it properly
    checked_eval_string("update = BMI.update");
    // Get function `update` and prepare arguments
    jl_function_t *func = jl_get_function(jl_main_module, "update");
    jl_value_t* arg1 = (jl_value_t*)jl_get_global(jl_main_module, jl_symbol("simulation"));
    jl_value_t* arg2 = jl_box_float64(dt);
    // Call function
    jl_call2(func, arg1, arg2);
    // Manual check if exception has been thrown
    check_exception();
    return 0;
}

ATHESIS_API int finalize()
{
     checked_eval_string("BMI.finalize(simulation)");
     return 0;
}

ATHESIS_API int get_start_time(double* start_time)
{
    jl_value_t *ret = checked_eval_string("BMI.get_start_time(simulation)");
    *start_time = jl_unbox_float64(ret);
    return 0;
}

ATHESIS_API int get_end_time(double* end_time)
{
    jl_value_t *ret = checked_eval_string("BMI.get_end_time(simulation)");
    *end_time = jl_unbox_float64(ret);
    return 0;
}

ATHESIS_API int get_current_time(double* current_time)
{
    jl_value_t *ret = checked_eval_string("BMI.get_current_time(simulation)");
    *current_time = jl_unbox_float64(ret);
    return 0;
}

ATHESIS_API int get_time_step(double* time_step)
{
    jl_value_t *ret = checked_eval_string("BMI.get_time_step(simulation)");
    *time_step = jl_unbox_float64(ret);
    return 0;
}

ATHESIS_API int set_var(char* var_name, void* var)
{    
    // TODO
    return 0;
}

ATHESIS_API int get_var(char* var_name, void* var)
{
    // TODO
    return 0;
}
