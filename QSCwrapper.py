from qsc import Qsc
from simsopt._core.graph_optimizable import Optimizable
import numpy as np

class QSCWrapper(Qsc, Optimizable):
    def __init__(self, *args, **kwargs):
        Qsc.__init__(self, *args, **kwargs)
        Optimizable.__init__(self, x0=Qsc.get_dofs(self),
                             external_dof_setter=Qsc.set_dofs,
                             names=self.names)

    def change_qsc_nfourier(self, ncoeffs):
        self.change_nfourier(ncoeffs)
        Optimizable.__init__(self, x0=Qsc.get_dofs(self),
                             external_dof_setter=Qsc.set_dofs,
                             names=self.names)

    def get_iota(self):
        return self.iota

    def get_max_elongation(self):
        return self.max_elongation

    def get_elongation(self):
        return self.elongation

    def get_B20QI_deviation(self):
        return self.B20QI_deviation

    def get_B2cQI_deviation(self):
        return self.B2cQI_deviation

    def get_B2sQI_deviation(self):
        return self.B2sQI_deviation

    def get_B20QI_deviation_max(self):
        return self.B20QI_deviation_max

    def get_B2cQI_deviation_max(self):
        return self.B2cQI_deviation_max

    def get_B2sQI_deviation_max(self):
        return self.B2sQI_deviation_max

    def get_X3c1(self):
        return self.X3c1

    def get_X3s1(self):
        return self.X3s1

    def get_Y3c1(self):
        return self.Y3c1

    def get_Y3s1(self):
        return self.Y3s1

    def get_delta(self):
        return self.delta

    def get_B0_well_depth(self):
        return self.B0_well_depth

    def get_k_second_order_SS(self):
        return self.k_second_order_SS

    def get_inv_L_grad_B(self):
        return self.inv_L_grad_B

    def get_B2c_array(self):
        return self.B2c_array

    def get_B2s_array(self):
        return self.B2s_array

    def get_X20(self):
        return self.X20

    def get_X2c(self):
        return self.X2c

    def get_X2s(self):
        return self.X2s

    def get_Y20(self):
        return self.Y20

    def get_Y2c(self):
        return self.Y2c

    def get_Y2s(self):
        return self.Y2s

    def get_Z20(self):
        return self.Z20

    def get_Z2c(self):
        return self.Z2c

    def get_Z2s(self):
        return self.Z2s

    def get_d_svals(self):
        return self.d_svals

    def get_max_d(self):
        return np.max(self.d)

    def get_d_curvature_d_varphi_at_0(self):
        return self.d_curvature_d_varphi_at_0

    def get_d_d_d_varphi_at_0(self):
        return self.d_d_d_varphi_at_0