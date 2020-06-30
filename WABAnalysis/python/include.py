"""Include this module"""
def library() :
    """Attach the name of WABAnalysis library to the process"""
    from LDMX.Framework.ldmxcfg import Process
    Process.addLibrary('/nfs/slac/g/ldmx/users/smidd/ldmx-sw/install/lib/libWABAnalysis.so')