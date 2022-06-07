import distutils.core
import Cython.Build
distutils.core.setup(
    ext_modules = Cython.Build.cythonize("TMP_1e_anomaly_ETo_step0_mod0_GMAO.pyx"))
