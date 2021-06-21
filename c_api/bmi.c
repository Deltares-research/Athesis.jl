// This file is a part of Julia. License is MIT: https://julialang.org/license

#include <julia.h>
#include <stdio.h>
#include <math.h>

#ifdef _OS_WINDOWS_
__declspec(dllexport) __cdecl
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

int main()
{
    jl_init();

    // {
    //     // Importing a Julia package

    //     checked_eval_string(
    //     "let dir = dirname(unsafe_string(Base.JLOptions().julia_bin))\n"
    //     // disable the package manager
    //     "    ENV[\"JULIA_PKGDIR\"] = joinpath(dir, \"disabled\")\n"
    //     // locate files relative to the "embedding" executable
    //     "    stdlib = filter(env -> startswith(Base.find_package(\"Distributed\"), env), Base.load_path())[end]\n"
    //     "    push!(empty!(LOAD_PATH), dir, stdlib)\n"
    //     "end"
    //     );
    //     checked_eval_string("import LocalModule");
    //     checked_eval_string("LocalModule.myapp()");
    // }

    // {
    //     // Main.include and Main.eval exist (#28825)
    //     checked_eval_string("include(\"include_and_eval.jl\")");
    //     checked_eval_string("f28825()");
    // }

    {
        // Test Athesis
        checked_eval_string("import Pkg");
        checked_eval_string("Pkg.activate(\"..\")");
        checked_eval_string("Pkg.instantiate()");
        checked_eval_string("using Athesis");
        checked_eval_string("include(\"/mnt/c/checkouts/Athesis/Athesis.jl/test/runtests.jl\")");

    }

    int ret = 0;
    jl_atexit_hook(ret);
    return ret;
}
