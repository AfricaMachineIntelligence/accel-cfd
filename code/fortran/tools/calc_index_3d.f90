      subroutine calc_index_3d(indx,i,j,k,m,n,l)
      
      integer, intent(in) :: i,j,k,m,n,l
      integer, intent(out) :: indx
      
      integer :: dm, dn, dl, 
      integer :: b_clc, b_lsl, b_sis
      integer :: plane1, plane2
      
      ! Create Abbriviation for dnterior lengths
      dm = m-2
      dn = n-2
      dl = l-2

      ! Create Abbrivation for blocks
      b_clc = 2*4+dm*5
      b_lsl = 2*5+dm*6
      b_sis = 2*6+dm*7
      plane1 = 2*b_clc+dn*b_lsl
      plane2 = 2*b_lsl+dn*b_sis

      ! Run Index Check and Calculate
      if (k == 1) then
          if (j == 1) then
              if (i == 1) then
                  indx = 1  
              else
                  indx = 5*(i-2)+4+2
              end if
          else if (j < n) then
              if (i == 1) then
                  indx = 2+(j-2)*b_lsl+b_clc
              else
                  indx = 6*(i-2)+8+(j-2)*b_lsl+b_clc
              end if
          else
              if (i == 1) then
                  indx = 2+dn*b_lsl+b_clc
              else
                  indx = 5*(i-2)+7+dn*b_lsl+b_clc
              end if
          end if
    
      else if (k < l) then
          if (j == 1) then
              if (i == 1) then
                  indx = 2+(k-2)*plane2+plane1  
              else
                  indx = 6*(i-2)+8+(k-2)*plane2+plane1
              end if
          else if (j < n) then
              if (i == 1) then
                  indx = 3+(j-2)*b_sis+b_lsl+(k-2)*plane2+plane1
              else
                  indx = 7*(i-2)+10+(j-2)*b_sis+b_lsl+(k-2)*plane2+plane1
              end if
          else
              if (i == 1) then
                  indx = 3+dn*b_sis+b_lsl+(k-2)*plane2+plane1
              else
                  indx = 6*(i-2)+9+dn*b_sis+b_lsl+(k-2)*plane2+plane1
              end if
          end if
    
      else
          if (j == 1) then
              if (i == 1) then
                  indx = 2+dl*plane2+plane1  
              else
                  indx = 5*(i-2)+7+dl*plane2+plane1
              end if
          else if (j < n) then
              if (i == 1) then
                  indx = 3+(j-2)*b_lsl+b_clc+dl*plane2+plane1
              else
                  indx = 6*(i-2)+9+(j-2)*b_lsl+b_clc+dl*plane2+plane1
              end if
          else
              if (i == 1) then
                  indx = 3+dn*b_lsl+b_clc+dl*plane2+plane1
              else
                  indx = 5*(i-2)+8+dn*b_lsl+b_clc+dl*plane2+plane1
              end if
          end if
      end if
      
      end subroutine calc_index_3d