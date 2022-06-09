
module prec
integer, parameter::hpi=1
integer, parameter::spi=2
integer, parameter::dpi=4
integer, parameter::mpi=4

integer, parameter::sp=selected_real_kind(p=6,r=37)
integer, parameter::dp=selected_real_kind(p=15,r=307)
integer, parameter::mp=selected_real_kind(p=15,r=307)

integer, parameter::spc=kind((1.0 , 1.0))
integer, parameter::dpc=kind((1.0_dp , 1.0_dp))
integer, parameter::mpc=kind((1.0_mp , 1.0_mp))

integer, parameter::lgt=kind(.true.)

end module prec
