#
# Wrapper script for Java Conda packages that ensures that the java runtime
# is invoked with the right options. Adapted from the bash script (http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128).

#
# Program Parameters
#
import os
import subprocess
import sys
import shutil
from os import access
from os import getenv
from os import X_OK
import os.path
from os import path

pkg_name = "transposon_annotation_tools_helitronscanner"
jar_file = "HelitronScanner.jar"
head_file = "head.lcvs"
tail_file = "tail.lcvs"

default_jvm_mem_opts = ['-Xms512m', '-Xmx1g']

# !!! End of parameter section. No user-serviceable code below this line !!!


def real_dirname(path):
    """Return the symlink-resolved, canonicalized directory-portion of path."""
    return os.path.dirname(os.path.realpath(path))


def java_executable():
    """Return the executable name of the Java interpreter."""
    java_home = getenv('JAVA_HOME')
    java_bin = os.path.join('bin', 'java')

    if java_home and access(os.path.join(java_home, java_bin), X_OK):
        return os.path.join(java_home, java_bin)
    else:
        return 'java'


def jvm_opts(argv):
    """Construct list of Java arguments based on our argument list.

    The argument list passed in argv must not include the script name.
    The return value is a 3-tuple lists of strings of the form:
      (memory_options, prop_options, passthrough_options)
    """
    mem_opts = []
    prop_opts = []
    pass_args = []
    exec_dir = None

    for arg in argv:
        if arg.startswith('-D'):
            prop_opts.append(arg)
        elif arg.startswith('-XX'):
            prop_opts.append(arg)
        elif arg.startswith('-Xm'):
            mem_opts.append(arg)
        elif arg.startswith('--exec_dir='):
            exec_dir = arg.split('=')[1].strip('"').strip("'")
            if not os.path.exists(exec_dir):
                shutil.copytree(real_dirname(sys.argv[0]), exec_dir, symlinks=False, ignore=None)
        else:
            pass_args.append(arg)

    # In the original shell script the test coded below read:
    #   if [ "$jvm_mem_opts" == "" ] && [ -z ${_JAVA_OPTIONS+x} ]
    # To reproduce the behaviour of the above shell code fragment
    # it is important to explictly check for equality with None
    # in the second condition, so a null envar value counts as True!

    if mem_opts == [] and getenv('_JAVA_OPTIONS') is None:
        mem_opts = default_jvm_mem_opts

    return (mem_opts, prop_opts, pass_args, exec_dir)


def def_temp_log_opts(args):
    """
    Establish default temporary and log folders.
    """
    TEMP  = os.getenv("TEMP")

    if TEMP is not None:
        if '-log' not in args:
            args.append('-log')
            args.append(TEMP+'/logs')

        if '-temp_folder' not in args :
            args.append('-temp_folder')
            args.append(TEMP)

    return args

def main():
    print("it was called")
    java = java_executable()
    """
    PeptideShaker updates files relative to the path of the jar file.
    In a multiuser setting, the option --exec_dir="exec_dir"
    can be used as the location for the peptide-shaker distribution.
    If the exec_dir dies not exist,
    we copy the jar file, lib, and resources to the exec_dir directory.
    """
    (mem_opts, prop_opts, pass_args, exec_dir) = jvm_opts(sys.argv[1:])
    pass_args = def_temp_log_opts(pass_args)
    jar_dir = exec_dir if exec_dir else real_dirname(sys.argv[0])

    if pass_args != [] and pass_args[0].startswith('eu'):
        jar_arg = '-cp'
    else:
        jar_arg = '-jar'
    
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in    
    #jar_path = os.path.join(jar_dir, jar_file)
    #jar_path = os.path.join(script_dir, '..', 'share', pkg_name, 'helitronscannerRES', 'jar', jar_file)    
    #hea_path = os.path.join(script_dir, '..', 'share', pkg_name, 'helitronscannerRES', 'lcvs', head_file)
    #tai_path = os.path.join(script_dir, '..', 'share', pkg_name, 'helitronscannerRES', 'lcvs', tail_file)
    jar_path = os.path.join(script_dir, pkg_name, 'helitronscannerRES', 'jar', jar_file)    
    hea_path = os.path.join(script_dir, pkg_name, 'helitronscannerRES', 'lcvs', head_file)
    tai_path = os.path.join(script_dir, pkg_name, 'helitronscannerRES', 'lcvs', tail_file)
    print("JarPath : ",jar_path)
    if(len(pass_args)==0):
        print("Please enter a command")
        print("For more help please type helitronscanner -help")
    elif(pass_args[0]=="-help"):
        print("HelitronScanner v1.0")
        print("Software can be run in three different modes...")
        print("  1) scanHead (scan DNA sequence for start points of transposons)")
        print("     helitronscanner scanHead -g dna.fasta -bs 0 -o headresult.txt")
        print("  2) scanTail (scan DNA sequence for end points of transposons)")
        print("     helitronscanner scanTail -g dna.fasta -bs 0 -o tailresult.txt")
        print("  3) pairends (group start and end points to identify helitron annotations)")
        print("     helitronscanner pairends -hs headresult.txt -ts tailresult.txt -o helitrons_result.txt")
        print("(default head.lcvs and tail.lcvs from author are used)")
    elif(pass_args[0]=="scanHead"):
        java_args = [java] + mem_opts + prop_opts + [jar_arg] + [jar_path] + ["scanHead", "-lf", hea_path] + pass_args[1:]
        sys.exit(subprocess.call(java_args))
    elif(pass_args[0]=="scanTail"):
        java_args = [java] + mem_opts + prop_opts + [jar_arg] + [jar_path] + ["scanTail", "-lf", tai_path] + pass_args[1:]
        sys.exit(subprocess.call(java_args))
    elif(pass_args[0]=="pairends"):
        java_args = [java] + mem_opts + prop_opts + [jar_arg] + [jar_path] + ["pairends"] + pass_args[1:]
        sys.exit(subprocess.call(java_args))
    else:
        print("HelitronScanner does not recognize command \"", pass_args[0],"\"...")
        print("For more help please type helitronscanner -help")
    #cmd = pass_args[0]
    #java_args = [java] + mem_opts + prop_opts + [jar_arg] + [jar_path] + pass_args
    #print("head exists...",path.exists(hea_path))
    #print("tail exists...",path.exists(tai_path))
    #sys.exit(subprocess.call(java_args))


if __name__ == '__main__':
    main()
