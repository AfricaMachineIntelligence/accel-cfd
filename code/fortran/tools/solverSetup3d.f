! This program is written in free-form

subroutine solverSetup3d (m, n, l, c, r, Ab, As, Aw, Ap, Ae, An, At, values, columns, rowIndex)

    integer, intent(in) :: m,n,l
    integer, intent(in) :: c, r
     
    real(dp), dimension(m, n, l), intent(in) :: Ab, As, Aw, Ap, Ae, An, At
    real(dp), dimension(c), intent(out) :: values, columns, 
    real(dp), dimension(r), intent(out) :: rowIndex

    integer, intent(in) :: vectorItr, valueItr, rowItr

    real(dp), dimension(7) :: valuesPlaceholder, columnsPlaceholder
    
    implicit none

!   Set initial values for loop iterators
    i = 1
    j = 1
    k = 1

!   Set initial values for index iterators
    vectorItr = 1
    valueItr = 1
    rowItr = 1

!   Initialize placeholder vectors
    valuesPlaceholder = 0
    columnsPlaceholder = 0

!   Initialize output vectors
    values = 0
    columns = 0
    rowIndex = 0

!   Start Fill-In Loop
    
    do k = 1, l

      do j = 1, n

        do i = 1,m

!         Check for node to the bottom
          if (k > 1) then

            columnsPlaceholder(vectorItr) = i + (j - 1) * m + (k - 2) * m * n
            valuesPlaceholder(vectorItr) = Ab(i,j,k)
            vectorItr = vectorItr + 1
          
          end if

!         Check for node to the south
          if (j > 1) then
          
            columnsPlaceholder(vectorItr) = i + (j - 2) * m + (k - 1) * m * n
            valuesPlaceholder(vectorItr) = As(i, j, k)
            vectorItr = vectorItr + 1
          
          end if

!         Check for node to the west
          if (i > 1) then

            columnsPlaceholder(vectorItr) = i - 1 + (j - 1) * m + (k - 1) * m * n
            valuesPlaceholder(vectorItr) = Aw(i, j, k)
            vectorItr = vectorItr + 1

          end if 

!         Update current node
          columnsPlaceholder(vectorItr) = i + (j - 1) * m + (k - 1) * m * n
          valuesPlaceholder(vectorItr) = Ap(i, j, k)
          vectorItr = vectorItr + 1

!         Check for node to the east
          if (i < m) then

            columnsPlaceholder(vectorItr) = i + 1 + (j - 1) * m + (k - 1) * m * n
            valuesPlaceholder(vectorItr) = Ae(i, j, k)
            vectorItr = vectorItr + 1

          end if

!         Check for node to the north
          if (j < n) then
          
            columnsPlaceholder(vectorItr) = i + j *m + (k - 1) * m * n
            valuesPlaceholder(vectorItr) = An(i, j, k)
            vectorItr = vectorItr + 1
          
          end if

!         Check for node to the top
          if (k < l) then

            columnsPlaceholder(vectorItr) = i + (j - 1) * m + k * m * n
            valuesPlaceholder(vectorItr) = At(i, j, k)
            vectorItr = vectorItr + 1

          end if

!         Copy placeholder values to output vectors
          values(valuesItr:valuesItr + (vectorItr - 2)) = valuesPlaceholder(1:vectorItr - 1)
          columns(valuesItr:valuesItr + (vectorItr - 2)) = columnsPlaceholder(1:vectorItr - 1)
          rowIndex(rowItr) = columnsPlaceholder(1)

!         Update index iterators
          rowItr = rowItr + 1
          valuesItr = valuesItr + vectorItr - 1
          vectorItr = 1

        end do

      end do

    end do

    return

end subroutine