!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine pulse(H0,Ht,time,count)
!***********************************************************************
! Purpose :  construct the time-dependent hamiltonian
!                        Ht = H0 + dipole*E_field(time)
!            argument 'count' is for storage of pulse
!            if count=-1 then current value of pulse is not stored
!***********************************************************************
  use rhodyn_data, only: zero, E_field, N_pulse, pulse_type,&
                         pulse_vec, pulse_vector, amp, taushift,&
                         sigma, power_shape, omega, phi,&
                         lu_pls, temp_vec, out_t, out_pulse,&
                         dipole_basis, t_local, omega_local,&
                         linear_chirp, flag_acorrection
  use mh5, only: mh5_put_dset
  use constants, only: pi, auToFs
  use definitions, only: iwp, wp
  implicit none
  complex(kind=wp), dimension(:,:), intent(in) :: H0
  complex(kind=wp), dimension(:,:), intent(out) :: Ht
  real(kind=wp), intent(in)       :: time
  integer(kind=iwp), intent(in)   :: count
  integer(kind=iwp) :: i

  E_field = zero

  do i=1,N_pulse
    t_local = time-taushift(i)
    omega_local = omega(i) + linear_chirp * t_local
    pulse_vec(:) = pulse_vector(i,:)
!
! sine^N pulse
! A\vec{e}\sin^n(\pi(t-t_0)/(2\sigma))\sin{(\Omega(t-t_0)+\varphi_0)}
! duration of the pulse equals 2\sigma
    if (pulse_type(1:3)=='SIN') then
      if (abs(t_local)<=sigma(i)) then
        E_field = E_field + amp(i)*pulse_vec &
        *sin( pi*t_local / (2.d0*sigma(i)) )**power_shape &
        *sin(omega_local * t_local + phi(i))
! correction due to vector potential derivative [Paramonov_JPCA_2012]:
! -A\vec{e}\pi/(2\sigma\Omega)
!     \sin(\pi(t-t_0)/(\sigma))\cos{(\Omega(t-t_0)+\varphi_0)}
! only for sine square shape form
        if (flag_acorrection.and.power_shape==2) then
          E_field = E_field - amp(i)*pulse_vec*pi &
          /(2*sigma(i)*omega_local) &
          *sin( pi*t_local / sigma(i) ) &
          *cos(omega_local * t_local + phi(i))
        endif
!      else
!        E_field=0d0
      endif
!
! cos^N pulse:
! A\vec{e}\cos^n(\pi(t-t_0)/(2\sigma))\sin{(\Omega(t-t_0)+\varphi_0)}
! duration of the pulse equals 2\sigma
    elseif (pulse_type(1:3)=='COS') then
      if (abs(t_local)<=sigma(i)) then
        E_field = E_field + amp(i)*pulse_vec &
        *cos( pi*t_local / (2.d0*sigma(i)) )**power_shape &
        *sin(omega_local * t_local + phi(i))
! correction is the same as in sine case, but with different choice of
! vector potential, here vector potential is given:
! -\vec{e}\frac{A}{\Omega} \cos^2(\pi(t-t_0)/(2\sigma))
!     \cos{(\Omega(t-t_0)+\varphi_0)}
! only for cosine square shape form
        if (flag_acorrection.and.power_shape==2) then
          E_field = E_field - amp(i)*pulse_vec*pi &
          /(2*sigma(i)*omega_local) &
          *sin( pi*t_local / sigma(i) ) &
          *cos(omega_local * t_local + phi(i))
        endif
!      else
!        E_field=0d0
      endif
!
! gaussian pulse
! A\vec{e}exp{-(t-t_0)^2/(2\sigma^2)}\sin{(\Omega(t-t_0)+\varphi_0)}
    elseif (pulse_type=='GAUSS') then
      E_field = E_field + amp(i)*pulse_vec &
      *exp(-t_local**2 / (2*sigma(i)**2)) &
      *sin(omega_local * t_local + phi(i))
! correction due to vector potential derivative:
! \vec{e}\frac{A(t-t_0)}{\sigma^2\Omega}
!     \exp{-(t-t_0)^2/(2\sigma)^2}\cos{(\Omega(t-t_0) + \varphi_0)}
      if (flag_acorrection) then
        E_field = E_field + amp(i)*pulse_vec*t_local &
        /(sigma(i)**2*omega_local) &
        *exp(-t_local**2 / (2*sigma(i)**2)) &
        *cos(omega_local * t_local + phi(i))
      endif
!
! monochromatic pulse
! A\vec{e}\sin{(\Omega(t-t_0)+\varphi_0)}
    elseif (pulse_type=='MONO') then
      E_field = amp(i)*pulse_vec &
      *sin(omega_local * (time-taushift(i)) + phi(i))
!
! explicitely polarized pulses
! think of more clever definition
    elseif (pulse_type=='MONO_R_CIRCLE') then
      E_field(1) = amp(i)*sin(omega(i)*time+phi(i))
      E_field(2) = Amp(1)*cos(omega(i)*time+phi(i))
      E_field(3) = 0d0
!
    elseif (pulse_type=='MONO_L_CIRCLE') then
      E_field(1) = amp(i)*cos(omega(i)*time+phi(i))
      E_field(2) = amp(i)*sin(omega(i)*time+phi(i))
      E_field(3) = 0d0
!
    elseif (pulse_type=='GAUSS_R_CIRCLE') then
      E_field(1) = amp(i)*sin(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(2) = amp(i)*cos(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(3)=0d0
!
    elseif (pulse_type=='GAUSS_L_CIRCLE')then
      E_field(1) = amp(i)*cos(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(2) = amp(i)*sin(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(3) = 0d0
    endif
  enddo
!
! saving pulse to files
  if (count>=0) then
    write(lu_pls,'(7(g25.15e3,2x))') time*auToFs, &
                  (dble(E_field(i)),aimag(E_field(i)), i=1,3)
    do i=1,3
      temp_vec(2*i-1) = dble (E_field(i))
      temp_vec(2*i)   = aimag(E_field(i))
    enddo
    call mh5_put_dset(out_pulse,temp_vec, [1,6],[count-1,0])
    call mh5_put_dset(out_t, [time*auToFs], [1], [count-1])
  endif
!
! update Hamiltonian H0 to Ht adding
! scalar product of E_field and dipole moment
  Ht = H0
  do i=1,3
    Ht(:,:) = Ht(:,:) + dipole_basis(:,:,i)*E_field(i)
  enddo
!
end
