# -*- coding: utf-8 -*-

#自动处理生成代码中的错误


import os
import sys

def file_replace(fname, s_0, s_1):
    print('replace %s with %s in %s' % (s_0, s_1, fname))
    if '\n' in s_0:
        f = open(fname,'r',encoding='utf-8')
        s = f.read()
        f.close()
        s = s.replace(s_0, s_1)
        f_new = open(fname,'w',encoding='utf-8')
        f_new.write(s)
        f_new.close()
        print("ok")
    else:
        fname_new = fname+'_'
        f = open(fname,'r',encoding='utf-8')
        f_new = open(fname_new,'w',encoding='utf-8')
        n = 0
        for line in f:
            if line.find(s_0) >= 0:
                n += 1
                line = line.replace(s_0, s_1)
            f_new.write(line)

        f.close()
        f_new.close()
        os.remove(fname)
        os.rename(fname_new, fname)
        print("found %d lines" % n)


def process(folder, pkgname):
    replace_list = [
    ['setup.py', 'cvxpy_codegen_solver', pkgname],
    ['setup.py', "'rt'", ''], # rt.lib
    ['setup.py', '\\', '/'], # bad slash
    ['codegenmodule.c', 'cvxpy_codegen_solver', pkgname],
    ['cvxpy_codegen_solver.py', 'cvxpy_codegen_solver', pkgname],
    ['codegen.h', '[0]', '[1]/*[zero]*/'], # zero len array
    ['codegen.h', 'vars_struct{\n}', 'vars_struct{int wtf;}'], # empty struct
    ['codegen.h', 'params_struct{}', 'params_struct{int wtf;}'], # empty struct
    ['ecos/src/ecos.c', '_set_output_format(_TWO_DIGIT_EXPONENT)', '0'], # old api
    ['ecos/include/glblopts.h', '#define NAN', '//#define _NAN'],
    ['ecos/include/glblopts.h', '#define INFINITY', '//#define _INFINITY'],
    ['param.c', '[0] = {};', '[1] = {0};'], # zero len array
    ]
    
    for replace_item in replace_list:
        fpath = os.path.join(folder, replace_item[0])
        if os.path.exists(fpath):
            file_replace(fpath, replace_item[1], replace_item[2])
    
    pyf = os.path.join(folder, 'cvxpy_codegen_solver.py')
    if os.path.exists(pyf):
        os.rename(pyf, os.path.join(folder, pkgname+'.py'))
    

if __name__ == '__main__':
    if len(sys.argv) > 2:
        folder = sys.argv[1].strip()
        pkgname = sys.argv[2].strip()
        process(folder, pkgname)
    else:
        print(sys.argv)
        print('missing par')