# Linux

1. Download and unpack julia in folder `~/Julia`
2. `cd c_api`
3. `~/Julia/julia-1.6.1/bin/julia ~/Julia/julia-1.6.1/share/julia/julia-config.jl --cflags --ldflags --ldlibs | xargs gcc -shared -o libathesis.so bmi.c`

# Windows
1. Download and install Julia
2. Open "MSYS2 MinGW 64-bit"
3. Go to "c_api"
4. `"$LOCALAPPDATA/Programs/Julia-1.6.1/bin/julia" "$LOCALAPPDATA/Programs/Julia-1.6.1/share/julia/julia-config.jl" --cflags --ldflags --ldlibs | xargs gcc -shared -o athesis.dll bmi.c`