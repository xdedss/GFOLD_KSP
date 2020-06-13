# GFOLD_KSP

连接krpc并调用求解器进行计算和自动降落的主要程序是demo3_gfold.py

可调参数在params.txt中

# codegen

运行GFOLD_codegen.py可以生成效率较高的c代码并作为python包安装到python环境里。

注意生成的代码可能有许多错误需要看情况修正，因此我写了postprocess.py来自动改正一些常见问题，但是不保证在所有环境下都能完美工作。

上述整个代码生成+改错+编译安装过程可以通过运行build.bat自动完成，如果一切正常，运行完之后python里就会多出两个包 gfold_solver_p3 和 gfold_solver_p4 

# 直接在python中求解

如果要直接cvxpy求解（耗时可能达到c的10倍以上）请把params.txt中direct改成True，然后应该就能直接运行了

# 测试

运行python ./GFOLD_run.py，可以用默认的测试数据+codegen的代码求解并显示结果（必须要先完成codegen的步骤）

加参数运行python ./GFOLD_run.py direct，则是同样的默认数据用直接求解的方法求解并显示结果



# License

此项目由[G-FOLD-Python](https://github.com/jonnyhyman/G-FOLD-Python)改进而成，遵循GPLv3
