!
! Module for handling particular organic matter (POM)
! To have other groups produce POM, do the following:
!  1) Define the masses of POM that each group produces in the vector mPOM 
!     This vector is already set to default as mPOM = m, which works for unicellular groups
!  2) Define the fluxes to POM in jPOM (note this is a rate 1/day)
!
! NOTE: THE DYNAMICS OF POM IN GENERALISTS AND COPEPODS NEEDS TO BE REVISITED
!
module POM
    use globals
    use spectrum
    use input
    implicit none
  
    private 
  
    !real(dp), parameter:: rhoCN = 5.68
    real(dp) :: remPOM ! = 0.07d0 ! remineralisation rate (1/day) (Serra-Pompei (2022)) @10 degrees
    real(dp) :: mPOMmin ! = 1e-9 ! Smallest POM mass
    real(dp) :: vel1, vel2
  
    type, extends(typeSpectrum) :: spectrumPOM
            
    contains
      procedure, pass :: initPOM
      procedure :: calcDerivativesPOM
      procedure :: printRates => printRatesPOM
    end type spectrumPOM
   
    public read_namelist_POM, initPOM, spectrumPOM, calcDerivativesPOM, printRatesPOM
  
  contains

  subroutine read_namelist_POM()
    integer :: file_unit,io_err

    namelist /input_POM / remPOM, mPOMmin, vel1, vel2

    call open_inputfile(file_unit, io_err)
        read(file_unit, nml=input_POM, iostat=io_err)
        call close_inputfile(file_unit, io_err)

  end subroutine read_namelist_POM

  subroutine initPOM(this, n, mMax)
    class(spectrumPOM):: this
    integer, intent(in):: n
    real(dp), intent(in):: mMax

    call read_namelist_POM()
    
    call this%initSpectrum(n, mPOMmin, mMax)

    this%velocity = vel1*this%m**vel2 ! Copepod fecal pellets from Serra-Pompei (2022)
    this%mort2 = 0.d0 ! No virulysis of POM
  end subroutine initPOM

  subroutine calcDerivativesPOM(this, u, dNdt, dDOCdt, dudt)
    class(spectrumPOM):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)
    real(dp):: uPOM(2)
    uPOM=u
    !write(*,*) " u(1) and velocity(1) is", u(1), this%velocity(1)
    !write(*,*) " u(2) and velocity(2) is", u(2), this%velocity(2)
    uPOM(1)=this%velocity(1)*remPOM*u(1)
    uPOM(2)=this%velocity(2)*remPOM*u(2)

    dudt = dudt - fTemp2*uPOM - this%mortpred*u
    dNdt = dNdt + fTemp2*sum(uPOM)/rhoCN
    dDOCdt = dDOCdt + fTemp2*sum(uPOM)
  end subroutine calcDerivativesPOM

  subroutine printRatesPOM(this)
    class(spectrumPOM), intent(in):: this
  
    write(*,*) "POM with ", this%n, " size classes:"
    call this%printRatesSpectrum()
  end subroutine printRatesPOM
end
