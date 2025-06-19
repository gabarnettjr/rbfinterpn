#!/usr/bin/perl

package Matrix;

# Methods and functions for creating and manipulating Matrix objects.
# Greg Barnett
# August 2023

use strict;
use warnings;

use constant PI => 3.141592653589793;

################################################################################

# TESTS

# /mu/bin/perl -e 'use lib "/home/gregorybarne/tmpCode"; use Matrix; Matrix::test_ALL;'

sub test_new
{
    my $A = Matrix::new([[1,2,3], [4,5,6], [7,8,9]]);
    print "test_new_\$A = \n";
    $A->disp;
}



sub test_rand
{
    my $A = Matrix::rand(2, 3);
    print "test_rand_\$A = \n";
    $A->disp;
}



sub test_alloc_set
{
    my $A = Matrix::alloc(4, 3);
    
    for (my $i = 0; $i < $A->numRows; $i++)
    {
        for (my $j = 0; $j < $A->numCols; $j++)
        {
            $A->set($i, $j, 1 + rand);
        }
    }
    
    print "test_alloc_set_\$A = \n";
    $A->disp;
}



sub test_zeros_ones_eye
{
    my $A = Matrix::zeros(4, 2);
    print "test_zeros_ones_eye_\$A = \n";
    $A->disp;
    
    my $B = Matrix::ones(3, 5);
    print "test_zeros_ones_eye_\$B = \n";
    $B->disp;
    
    my $c = Matrix::eye(5);
    print "test_zeros_ones_eye_\$C = \n";
    $c->disp;
}



sub test_linspace
{
    my $x = Matrix::linspace(-1, 3, 5);
    print "test_linspace_\$x = ";
    $x->disp;
}



sub test_round
{
    my $A = Matrix::rand(3, 7);
    print "test_round_\$A = \n";
    $A->disp;
    
    my $B = $A->round(2);
    print "test_round_\$B = \n";
    $B->disp;
}



sub test_rows_cols
{
    my $A = Matrix::rand(5, 3);
    print "test_rows_cols_\$A = \n";
    $A->disp;
    
    my $B = $A->rows([1,4]);
    print "test_rows_cols_\$B = \n";
    $B->disp;
    
    my $C = $A->cols([0,2]);
    print "test_rows_cols_\$C = \n";
    $C->disp;
}



sub test_meshgrid
{
    my $x = Matrix::linspace(-3, 3, 7);
    my $y = Matrix::linspace(1, 5, 5);
    print "test_meshgrid_\$x = \n";
    $x->disp;
    print "test_meshgrid_\$y = \n";
    $y->disp;
    my ($xx, $yy) = Matrix::meshgrid($x, $y);
    print "test_meshgrid_\$xx = \n";
    $xx->disp;
    print "test_meshgrid_\$yy = \n";
    $yy->disp;
}



sub test_vstack
{
    my $a = Matrix::rand(3, 2);
    my $b = Matrix::ones(4, 2);
    my $c = $a->vstack($b);
    print "test_vstack_\$c = \n";
    $c->disp;
}



sub test_hstack
{
    my $a = Matrix::rand(2, 3);
    my $b = Matrix::ones(2, 4);
    my $c = $a->hstack($b);
    print "test_hstack_\$c = \n";
    $c->disp;
}



sub test_flatten
{
    my $a = Matrix::rand(3, 5);
    my $b = $a->flatten;
    print "test_flatten_\$b = \n";
    $b->disp;
}



sub test_setRows
{
    my $A = Matrix::zeros(5, 5);
    print "test_setRows_\$A = \n";
    $A->disp;
    $A->setRows([0,4], Matrix::ones(2, 5));
    print "test_setRows_\$A = \n";
    $A->disp;
}



sub test_setCols
{
    my $A = Matrix::zeros(5, 5);
    print "test_setCols_\$A = \n";
    $A->disp;
    $A->setCols([0,4], Matrix::ones(5, 2));
    print "test_setCols_\$A = \n";
    $A->disp;
}



sub test_spliceRows
{
    my $A = Matrix::rand(5, 3);
    print "test_spliceRows_\$A = \n";
    $A->disp;
    $A = $A->spliceRows(1, 3);
    print "test_spliceRows_\$A = \n";
    $A->disp;
}



# sub test_spliceCols
# {
    # my $A = Matrix::rand(3, 7);
    # print "test_spliceCols_\$A = \n";
    # $A->disp;
    # $A = $A->spliceCols(1, 3);
    # print "test_spliceCols_\$A = \n";
    # $A->disp;
# }



sub test_norm
{
    my $A = Matrix::ones(1, 4)->dot(-1)->plus(Matrix::rand(1, 4)->dot(2));
    print "test_norm_\$A = \n";
    $A->disp;
    print "inf_norm(\$A) = " . $A->norm("inf") . "\n";
}



sub test_solve
{
    my $A = Matrix::eye(5);
    my $b = Matrix::ones(5, 1);
    my $x = $A->solve($b);
    
    print "test_solve_\$A = \n";
    $A->disp;
    print "test_solve_\$b = \n";
    $b->disp;
    print "test_solve_\$x = \n";
    $x->disp;
}



sub test_leastSquares
{
    my $A = Matrix::new([[0.1888460, 0.0394834, 0.0184743, 0.8008640, 0.818357]
                       , [0.4792610, 0.8637790, 0.6874270, 0.9958450, 0.272958]
                       , [0.0704303, 0.5977300, 0.0605322, 0.5260480, 0.771405]
                       , [0.3555720, 0.5175610, 0.3917400, 0.8764970, 0.629965]
                       , [0.0152446, 0.4146580, 0.9055900, 0.2382900, 0.232571]
                       , [0.5353680, 0.1206550, 0.0684741, 0.7197770, 0.718718]
                       , [0.1702700, 0.6308350, 0.9088270, 0.4175410, 0.269781]
                       , [0.8023090, 0.7830920, 0.5452130, 0.6701250, 0.505405]
                       , [0.8550240, 0.0936718, 0.7088540, 0.0640175, 0.209113]
                       , [0.8814950, 0.9379720, 0.4306120, 0.0812095, 0.843693]]);
    
    my $b = Matrix::new([[0.71376695596717880]
                       , [0.47569193631426410]
                       , [0.10628212554412164]
                       , [0.05427698291996352]
                       , [0.67411665587272970]
                       , [0.09467686960996125]
                       , [0.16625878478066125]
                       , [0.28295185053589456]
                       , [0.34237011509123440]
                       , [0.64927664119453020]]);
    
    my $x = $A->leastSquares($b);
    
    my $juliaAnswer = Matrix::new([[-0.03328356882014103]
                                 , [-0.14599611906283225]
                                 , [ 0.48350091879301293]
                                 , [-0.12348853474086356]
                                 , [ 0.53860316304772900]]);
    
    print "test_leastSquares_\$A = \n";
    $A->disp;
    print "test_leastSquares_\$b = \n";
    $b->disp;
    print "test_leastSquares_\$x = \n";
    $x->disp;
    print "test_leastSquares_\$juliaAnswer = \n";
    $juliaAnswer->disp;
}



sub test_sum
{
    my $A = Matrix::ones(5, 3);
    my $sum = $A->sum;
    
    print "test_sum_\$A = \n";
    $A->disp;
    print "\$sum = $sum\n\n";
}



sub test_avg
{
    my $A = Matrix::rand(10, 10);
    my $avg = $A->avg;
    
    print "test_avg_\$A = \n";
    $A->disp;
    print "\$avg = $avg\n\n";
}



sub test_ALL
{
    test_new;
    test_rand;
    test_alloc_set;
    test_zeros_ones_eye;
    test_linspace;
    test_round;
    test_rows_cols;
    test_meshgrid;
    test_vstack;
    test_hstack;
    test_flatten;
    test_setRows;
    test_setCols;
    test_spliceRows;
    # test_spliceCols;
    test_norm;
    test_solve;
    test_leastSquares;
    test_sum;
    test_avg;
}

################################################################################

# CONSTRUCTORS

sub new
{
    # Create new matrix from a 2D array of values.
    my %self;
    my $self = \%self;
    bless $self;
    
    if ((scalar @_ == 1 || scalar @_ == 2) && ref $_[0] && ref @{$_[0]}[0])
    {
        $$self{'type'}    = 'Matrix';
        $$self{'items'}   = shift;
        $$self{'numRows'} = scalar @{$$self{'items'}};
        $$self{'numCols'} = scalar @{@{$$self{'items'}}[0]};
        if (scalar @_)
        {
            $$self{'tol'} = shift;
        }
        else
        {
            $$self{'tol'} = 1e-8;                                      # DEFAULT
        }
        
        foreach my $row (@{$$self{'items'}})
        {
            if (scalar @{$row} != $$self{'numCols'})
            {
                print "\nEach row should have the same number of elements.\n";  die;
            }
        }
    }
    else
    {
        print "\nBad input for new matrix.  Give ref to array of refs to arrays.\n";  die;
    }
    
    return $self;
}



sub alloc
{
    # Allocate memory to fill in the items of a new matrix.
    # If all items will be filled, then this is faster than using zeros().
    if (scalar @_ != 2)
    {
        print "\nMatrix::alloc() requires two inputs: numRows and numCols.\n";  die;
    }
    my ($numRows, $numCols) = @_;
    
    my @items;
    $#items = ($numRows - 1);
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        my @tmp;
        $#tmp = ($numCols - 1);
        $items[$i] = \@tmp;
    }

    return Matrix::new(\@items);
}



sub zeros
{
    # New matrix of zeros.
    if (scalar @_ != 2)
    {
        print "\nMatrix::zeros() requires two inputs: numRows and numCols.\n";  die;
    }
    my ($numRows, $numCols) = @_;
    
    my @items = ();
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        my @tmp = (0) x $numCols;
        push(@items, \@tmp);
    }
    
    return Matrix::new(\@items);
}



sub ones
{
    # New matrix of ones.
    if (scalar @_ != 2)
    {
        print "\nMatrix::ones() requires two inputs: numRows and numCols.\n";  die;
    }
    my ($numRows, $numCols) = @_;
    
    my @items = ();
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        my @tmp = (1) x $numCols;
        push(@items, \@tmp);
    }
    
    return Matrix::new(\@items);
}



sub rand
{
    # New matrix whose values are randomly chosen from between 0 and 1.
    if (scalar @_ != 2)
    {
        print "\nMatrix::rand() requires two inputs: numRows and numCols.\n";  die;
    }
    my ($numRows, $numCols) = @_;
    
    my @items = ();
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        my @tmp = ();
        for (my $j = 0; $j < $numCols; $j++)
        {
            push(@tmp, rand);
        }
        push(@items, \@tmp);
    }
    
    return Matrix::new(\@items);
}



sub eye
{
    # New Identity matrix.
    if (scalar @_ != 1)
    {
        print "\nMatrix::eye() requires one input: numRows (numCols always equals numRows).\n";  die;
    }
    my ($numRows) = @_;

    my $numCols = $numRows;
    my @items = ();
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        my @tmp = ();
        for (my $j = 0; $j < $numCols; $j++)
        {
            if ($j == $i) 
            {
                push(@tmp, 1);
            }
            else
            {
                push(@tmp, 0);
            }
        }
        push(@items, \@tmp);
    }
    
    return Matrix::new(\@items);
}



sub linspace
{
    # New 1D matrix of equally-spaced values with known start and finish.
    if (scalar @_ != 3)
    {
        print "\nMatrix::linspace() requires three inputs: start, finish, and number of points.\n";  die;
    }
    my ($a, $b, $numCols) = @_;
    my $tmp = Matrix::new([[$numCols]]);
    if ($numCols < 1 || ! $tmp->round(0)->equals($tmp))
    {
        print "\nLast input must be a positive whole number greater than 0 (number of points in array).  Your input is \"$numCols\".\n";  die;
    }
    
    my $dx = 0;
    $dx = ($b - $a) / ($numCols - 1) if $numCols != 1;
    my $item = $a;
    my @items = ();
    for (my $i = 0; $i < $numCols; $i++)
    {
        push @items, $item;
        $item += $dx;
    }
    my $numRows = 1;
    
    return Matrix::new([\@items]);
}



sub round
{
    # Round the entries in a matrix to the specified number of places.
    my $self = shift;
    if (scalar @_ != 1)
    {
        print "\nOne input is required (number of places to round to).\n";  die;
    }
    my $places = shift;
    
    my $out = Matrix::alloc($self->numRows, $self->numCols);
    
    for (my $i = 0; $i < $out->numRows; $i++) {
        for (my $j = 0; $j < $out->numCols; $j++)
        {
            $out->set($i, $j, int (($self->item($i, $j) * 10**$places + .5)) / 10**$places) if $self->item($i, $j) >= 0;
            $out->set($i, $j, int (($self->item($i, $j) * 10**$places - .5)) / 10**$places) if $self->item($i, $j) <= 0;
        }
    }
    
    return $out;
}



sub equals
{
    my $self = shift;
    if (scalar @_ != 1)
    {
        print "\nExactly one input is required (other matrix).\n";  die;
    }
    my $other = shift;
    if ($self->numRows != $other->numRows || $self->numCols != $other->numCols)
    {
        print "\nTwo matrices must be same size to check if they are equal.\n";  die;
    }
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            my $absDiff = abs ($self->item($i, $j) - $other->item($i, $j));
            
            if ($absDiff > $self->tol || $absDiff > $other->tol)
            {
                return 0;
            }
        }
    }
    
    return 1;
}

################################################################################

# GETTERS

sub type
{
    my $self = shift;
    return $$self{'type'};
}



sub items
{
    my $self = shift;
    return $$self{'items'};
}



sub numRows
{
    my $self = shift;
    return scalar @{$$self{'items'}};
}



sub numCols
{
    my $self = shift;
    return scalar @{@{$$self{'items'}}[0]};
}



sub tol
{
    my $self = shift;
    return $$self{'tol'};
}



sub rows
{
    # Get the rows described in the ref to array $ind.
    my $self = shift;
    my $ind = shift;
    
    foreach my $i (@{$ind})
    {
        if ($i >= $self->numRows)
        {
            print "\nOne or more requested row index is out of range.\n";  die;
        }
    }
    
    my $out = Matrix::alloc(scalar @{$ind}, $self->numCols);
    
    for (my $i = 0; $i < $out->numRows;  $i++)
    {
        for (my $j = 0; $j < $out->numCols; $j++)
        {
            $out->set($i, $j, $self->item(@{$ind}[$i], $j));
        }
    }
    
    return $out;
}



sub cols
{
    # Get the columns described in the ref to array $ind.
    my $self = shift;
    my $ind = shift;
    
    foreach my $j (@{$ind})
    {
        if ($j >= $self->numCols)
        {
            print "\nOne or more requested column index is out of range.\n";  die;
        }
    }
    
    my $out = Matrix::alloc($self->numRows, scalar @{$ind});
    
    for (my $i = 0; $i < $out->numRows;  $i++)
    {
        for (my $j = 0; $j < $out->numCols; $j++)
        {
            $out->set($i, $j, $self->item($i, @{$ind}[$j]));
        }
    }
    
    return $out;
}



sub spliceRows
{
    my $self = shift;
    my $ind = shift;
    my $howMany = shift;
    
    my $out = $self->copy;
    splice @{$out->items}, $ind, $howMany;
    return $out;
}



# sub spliceCols
# {
    # # DANGER: This method is slow.  Transpose the matrix and use spliceRows.
    # my $self = shift;
    # my $ind = shift;
    # my $howMany = shift;
    
    # my $out = Matrix::alloc($self->numRows, $self->numCols - $howMany);
    # my $indices;
    
    # if ($ind > 0)
    # {
        # $indices = Matrix::linspace(0, $ind - 1, $ind);
        # $out->setCols(@{$indices->items}[0], $self->cols(@{$indices->items}[0]));
    # }
        
    # if ($ind + ($howMany - 1) < $self->numCols - 1)
    # {
        # $indices = Matrix::linspace($ind, $out->numCols - 1, $self->numCols - $ind - $howMany);
        # $out->setCols(@{$indices->items}[0], $self->cols(@{$indices->plus($howMany)->items}[0]));
    # }
    
    # return $out;
# }

################################################################################

sub len
{
    my $self = shift;
    return $self->numRows if $self->numCols == 1;
    return $self->numCols if $self->numRows == 1;
    return $self->numRows * $self->numCols;
}



sub meshgrid
{
    my $x = shift;
    my $y = shift;
    
    my $xx = Matrix::alloc($y->numCols, $x->numCols);
    my $yy = Matrix::alloc($y->numCols, $x->numCols);
    
    for (my $i = 0; $i < $y->numCols; $i++)
    {
        for (my $j = 0; $j < $x->numCols; $j++)
        {
            $xx->set($i, $j, $x->item($j));
            $yy->set($i, $j, $y->item($i));
        }
    }
    
    return ($xx, $yy);
}



sub vstack
{
    # Stack two matrices vertically.
    my $self = shift;
    my $other = shift;
    
    if ($self->numCols != $other->numCols)
    {
        print "\nDimension mismatch.\n";  die;
    }
    
    my $out = Matrix::alloc($self->numRows + $other->numRows, $other->numCols);
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $out->set($i, $j, $self->item($i, $j));
        }
    }
    
    for (my $i = $self->numRows; $i < $self->numRows + $other->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $out->set($i, $j, $other->item($i - $self->numRows, $j));
        }
    }
    
    return $out;
}



sub hstack
{
    # Stack to matrices horizontally.
    my $self = shift;
    my $other = shift;
    
    return $self->transpose->vstack($other->transpose)->transpose;
}



sub flatten
{
    # Turn the matrix into one long row by placing rows next to each other.
    my $self = shift;
    
    my $out = Matrix::alloc(1, $self->numRows * $self->numCols);
    my $k = 0;
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $out->set($k, $self->item($i, $j));
            $k++;
        }
    }
    
    return $out;
}



sub reshape
{
    # Reshape a row-vector to be the desired size.
    
    my $self = shift;

    if ((scalar @_) != 2)
    {
        print "\nExactly two inputs (numRows and numCols) are required.\n";  die;
    }

    my $numRows = shift;
    my $numCols = shift;

    if ($self->numRows * $self->numCols != $numRows * $numCols)
    {
        print "\nTo reshape, the total number of elements must equal numRows * numCols.\n";  die;
    }

    my $out = Matrix::alloc($numRows, $numCols);
    my $k = 0;
    my $flat = $self->flatten;

    for (my $i = 0; $i < $numRows; $i++)
    {
        for (my $j = 0; $j < $numCols; $j++)
        {
            $out->set($i, $j, $flat->item($k));
            $k++;
        }
    }

    return $out;
}



sub disp
{
    # Display the contents of a matrix in a readable format.
    my $self = shift;
    
    foreach my $row (@{$self->items})
    {
        foreach my $item (@{$row})
        {
            printf "%8.5f ", $item;
            # printf "%13.9f ", $item;
        }
        print "\n";
    }
    
    print "\n";
}



sub item
{
    # Get a single element of a matrix at row $i and column $j.
    my $self = shift;
    my ($i, $j);
    
    if (scalar @_ == 2)
    {
        # Matrix with rows and columns.
        $i = shift;
        $j = shift;
        
        if ($i >= $self->numRows || $j >= $self->numCols)
        {
            print "\nRow or column index out of range.\n";  die;
        }
        
        return @{@{$self->items}[$i]}[$j];
    }
    elsif (scalar @_ == 1 && $self->numRows == 1)
    {
        # Row vector.
        $j = shift;
        
        if ($j >= $self->numCols)
        {
            print "\nColumn index out of range.\n";  die;
        }
        
        return @{@{$self->items}[0]}[$j];
    }
    elsif (scalar @_ == 1 && $self->numCols == 1)
    {
        # Column vector.
        $i = shift;
        
        if ($i >= $self->numRows)
        {
            print "\nRow index out of range.\n";  die;
        }
        
        return @{@{$self->items}[$i]}[0];
    }
}



sub set
{
    # Set a single element of a matrix equal to a specified value.
    my $self = shift;    
    my ($i, $j, $val);
    
    # Undefine min/max info since something in the matrix is changing.
    if (defined $$self{'min'})
    {
        undef $$self{'min'};
        undef $$self{'iMin'};
        undef $$self{'jMin'};
        undef $$self{'max'};
        undef $$self{'iMax'};
        undef $$self{'jMax'};
    }
    
    if (scalar @_ == 3)
    {
        # Matrix with multiple rows and columns.
        $i = shift;
        $j = shift;
        
        if ($i >= $self->numRows || $j > $self->numCols)
        {
            print "\nRow or column index out of range.\n";  die;
        }
        
        $val = shift;
        @{@{$self->items}[$i]}[$j] = $val;
    }
    elsif (scalar @_ == 2 && $self->numRows == 1)
    {
        # Row vector.
        $j = shift;
        
        if ($j > $self->numCols)
        {
            print "\nColumn index out of range.\n";  die;
        }
        
        $val = shift;
        @{@{$self->items}[0]}[$j] = $val;
    }
    elsif (scalar @_ == 2 && $self->numCols == 1)
    {
        # Column vector.
        $i = shift;
        
        if ($i >= $self->numRows)
        {
            print "\nRow index out of range.\n";  die;
        }
        
        $val = shift;
        @{@{$self->items}[$i]}[0] = $val;
    }
    else
    {
        print "\nInputs not understood.\n";  die;
    }
}



sub setRows
{
    # Set certain rows of an existing matrix equal to new values.
    my $self = shift;
    my $ind = shift;
    my $vals = shift;
    
    foreach my $i (@{$ind})
    {
        if ($i >= $self->numRows)
        {
            print "\nAttempt to access row index that is out of range.\n";  die;
        }
    }
    
    for (my $i = 0; $i < scalar @{$ind}; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $self->set(@{$ind}[$i], $j, $vals->item($i, $j));
        }
    }
}



sub setCols
{
    # Set certain columns of an existing matrix equal to new values.
    my $self = shift;
    my $ind = shift;
    my $vals = shift;
    
    foreach my $j (@{$ind})
    {
        if ($j >= $self->numCols)
        {
            print "\nAttempt to access column index that is out of range.\n";  die;
        }
    }
    
    for (my $j = 0; $j < scalar @{$ind}; $j++)
    {
        for (my $i = 0; $i < $self->numRows; $i++)
        {
            $self->set($i, @{$ind}[$j], $vals->item($i, $j));
        }
    }
}



# sub containsCol
# {
    # my $self = shift;
    # my $other = shift;
    
    # my $tol = $self->tol;
    # $tol = $other->tol if $other->tol < $self->tol;
    
    # for (my $j = 0; $j < 
# }



sub copy
{
    my $self = shift;    
    my $out = Matrix::alloc($self->numRows, $self->numCols);
    
    for (my $i = 0; $i < $out->numRows; $i++)
    {
        for (my $j = 0; $j < $out->numCols; $j++)
        {
            $out->set($i, $j, $self->item($i, $j));
        }
    }
    
    return $out;
}

################################################################################

# OPERATIONS

sub plus
{
    # Add two matrices or add a scalar to a matrix.
    my $self = shift;
    my $other = shift;
    
    my $numRows = $self->numRows;
    my $numCols = $self->numCols;
    
    if ((ref $other) && ($numRows != $other->numRows || $numCols != $other->numCols)) {
        print "\nMatrices must be the same size to add them together.\n";  die;
    }
    
    my $sum = Matrix::alloc($numRows, $numCols);
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        for (my $j = 0; $j < $numCols; $j++)
        {
            if (ref $other)
            {
                $sum->set($i, $j, $self->item($i, $j) + $other->item($i, $j));
            }
            else
            {
                $sum->set($i, $j, $self->item($i, $j) + $other);
            }
        }
    }
    
    return $sum;
}



sub minus
{
    # Subtract two matrices or subtract a scalar from a matrix.
    my $self = shift;
    my $other = shift;

    if (ref $other)
    {
        return $self->plus($other->dot(-1));
    }
    else
    {
        return $self->plus(-1 * $other);
    }
}



sub dotProduct
{
    # Return the dot product of two 1D matrices of the same size.
    my $self = shift;
    my $other = shift;

    if (! ref $other)
    {
        print "\nInput must be a 1D matrix.\n";  die;
    }
    
    my $numRows = $self->numRows;
    my $numCols = $self->numCols;

    my $dotProduct = 0;
    
    if ($numRows != $other->numRows || $numCols != $other->numCols)
    {
        print "\nMatrices must be the same size to dot them.\n";  die;
    }
    elsif ($numRows != 1 && $numCols != 1)
    {
        print "\nThis function is only implemented for 1D matrices.\n";  die;
    }
    elsif ($numRows == 1)
    {
        for (my $j = 0; $j < $numCols; $j++)
        {
            $dotProduct += $self->item(0, $j) * $other->item(0, $j);
        }
    }
    elsif ($numCols == 1)
    {
        for (my $i = 0; $i < $numRows; $i++)
        {
            $dotProduct += $self->item($i, 0) * $other->item($i, 0);
        }
    }
    
    return $dotProduct;
}



sub transpose
{
    # Transpose a matrix so rows become columns and vice versa.
    my $self = shift;
    
    my $out = Matrix::alloc($self->numCols, $self->numRows);
    
    for (my $i = 0; $i < $out->numRows; $i++)
    {
        for (my $j = 0; $j < $out->numCols; $j++)
        {
            $out->set($i, $j, $self->item($j, $i));
        }
    }
    
    return $out;
}



sub dot
{
    # Multiply a matrix by a scalar, or by another matrix.
    my $self = shift;
    my $other = shift;
    
    my $prod;
    
    if (ref $other)
    {
        if ($self->numCols != $other->numRows)
        {
            print "\n\$numCols of first matrix (" . $self->numCols . ") must equal \$numRows of second (" . $other->numRows . ").\n";
            print "\$first = \n";
            $self->disp;
            print "\n";
            print "\$second = \n";
            $other->disp;
            die;
        }
        
        $prod = Matrix::alloc($self->numRows, $other->numCols);
        $other = $other->transpose;
        
        for (my $i = 0; $i < $prod->numRows; $i++) {
            for (my $j = 0; $j < $prod->numCols; $j++) {
                $prod->set($i, $j, $self->rows([$i])->dotProduct($other->rows([$j])));
            }
        }
    }
    else
    {
        $prod = Matrix::alloc($self->numRows, $self->numCols);
        
        for (my $i = 0; $i < $self->numRows; $i++)
        {
            for (my $j = 0; $j < $self->numCols; $j++)
            {
                $prod->set($i, $j, $self->item($i, $j) * $other);
            }
        }
    }
    
    return $prod;
}



sub dotTimes
{
    # Multiply two matrices element-wise.
    my $self = shift;
    my $other = shift;
    
    my $numRows = $self->numRows;
    my $numCols = $self->numCols;
    
    if (! ref $other)
    {
        print "\nUse dot() to multiply a scalar by a matrix.\n";  die;
    }
    elsif ($numRows != $other->numRows || $numCols != $other->numCols)
    {
        print "\nMatrices must be the same size to (dot) multiply them together.\n";  die;
    }
    
    my $prod = Matrix::alloc($numRows, $numCols);
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        for (my $j = 0; $j < $numCols; $j++)
        {
            $prod->set($i, $j, $self->item($i, $j) * $other->item($i, $j));
        }
    }
    
    return $prod;
}



sub dotDiv
{
    # Divide two matrices element-wise.
    my $self = shift;
    my $other = shift;
    
    my $numRows = $self->numRows;
    my $numCols = $self->numCols;
    
    if (! ref $other)
    {
        print "\nInput must be a matrix to divide by element-wise.\n";  die;
    }
    elsif ($numRows != $other->numRows || $numCols != $other->numCols)
    {
        print "\nMatrices must be the same size to (dot) divide them.\n";  die;
    }
    
    my $quot = Matrix::alloc($numRows, $numCols);
    
    for (my $i = 0; $i < $numRows; $i++)
    {
        for (my $j = 0; $j < $numCols; $j++)
        {
            $quot->set($i, $j, $self->item($i, $j) / $other->item($i, $j));
        }
    }
    
    return $quot;
}



sub pow
{
    # Raise a matrix to a power element-wise.
    my $self = shift;
    my $pow = shift;
    
    my $out = $self->copy();
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $out->set($i, $j, $self->item($i, $j) ** $pow);
        }
    }
    
    return $out;
}



sub sum
{
    # Return the sum of all elements of a matrix.
    my $self = shift;
    
    my $sum = 0;
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $sum += $self->item($i, $j);
        }
    }
    
    return $sum;
}



sub avg
{
    # Return the average of all elements of a matrix.
    my $self = shift;
    
    return $self->sum / $self->len;
}



sub std
{
    # Return the sample standard deviation of all values in a matrix.
    my $self = shift;
    
    my $out = 0;
    my $avg = $self->avg;
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $out += ($self->item($i, $j) - $avg) ** 2;
        }
    }
    
    $out = $out / ($self->len - 1);
    $out = sqrt $out;
    return $out;
}



sub mode
{
    my $self = shift;
    my @numSame = ();
    my @modes = ();
    my $flatSelf = $self->flatten;
    
    foreach my $item (@{@{$flatSelf->items}[0]})
    {
        if (! @numSame && ! @modes)
        {
            push @numSame, 1;
            push @modes, $item;
        }
        else
        {
            my $foundSame = 0;
            
            for (my $i = 0; $i < scalar @modes; $i++)
            {
                if (abs ($item - $modes[$i]) < $self->tol)
                {
                    $foundSame = 1;
                    $numSame[$i]++;
                }
            }
            
            if (! $foundSame)
            {
                push @modes, $item;
                push @numSame, 1;
            }
        }
        
        # # DEBUG
        # print "modes:   ";
        # foreach my $mode (@modes)
        # {
            # print $mode . " ";
        # }
        # print "\n";
        # print "numSame: ";
        # foreach my $same (@numSame)
        # {
            # print $same . " ";
        # }
        # print "\n\n";
    }
    
    # Return the mode.
    my $numSame = Matrix::new([\@numSame]);
    return ($numSame->max, $modes[$numSame->jMax]);
}



sub norm
{
    # Calculate the norm of a matrix (default is 2-norm).
    my $self = shift;
    if ($self->numRows > 1 && $self->numCols > 1)
    {
        print "\nOnly 1D matrices are supported right now.\n";  die;
    }
    my $p = shift;

    $p = 2 if ! $p;
    my $ell;
    
    if ($self->numCols == 1)
    {
        $ell = $self->numRows;
    }
    else
    {
        $ell = $self->numCols;
    }

    my $norm = 0;

    if ($p =~ /^inf(inity)?$/i)
    {
        $norm = abs $self->item(0);
        
        for (my $j = 1; $j < $ell; $j++)
        {
            my $val = abs $self->item($j);
            $norm = $val if $val > $norm;
        }
    }
    else
    {
        for (my $j = 0; $j < $ell; $j++)
        {
            $norm += (abs $self->item($j)) ** $p;
        }
        
        $norm = $norm ** (1 / $p);
    }
    
    return $norm;
}



sub cos
{
    my $self = shift;
    
    my $out = Matrix::alloc($self->numRows, $self->numCols);
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $out->set($i, $j, CORE::cos $self->item($i, $j));
        }
    }
    
    return $out;
}



sub sin
{
    my $self = shift;
    
    my $out = Matrix::alloc($self->numRows, $self->numCols);
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            $out->set($i, $j, CORE::sin $self->item($i, $j));
        }
    }
    
    return $out;
}



sub tan
{
    my $self = shift;
    return $self->sin->dotDiv($self->cos);
}



sub wiggle
{
    # Move the coordinates in the matrix in a random direction by a random
    # amount, where the max possible distance moved is controlled by $drMax.
    my $self = shift;
    if ((scalar @_) !~ /^(1|5)$/)
    {
        print "\nOne input (max wiggle amount) is required, or five inputs (also specify bounding box for rectangular coords).\n";  die;
    }
    if ($self->numRows != 2)
    {
        print "\nInput should be a coordinate matrix (two rows).\n";  die;
    }
    my $drMax = shift;
    
    # Optional inputs for coordinates in a rectangular domain with known corners.
    my $xll = shift;
    my $yll = shift;
    my $xur = shift;
    my $yur = shift;
    
    my $dr  = Matrix::rand(1, $self->numCols)->dot($drMax);
    my $dth = Matrix::rand(1, $self->numCols)->dot(2 * Matrix::PI);
    
    my $dx = $dr->dotTimes($dth->cos);
    my $dy = $dr->dotTimes($dth->sin);
    
    if (defined $xll && defined $yll && defined $xur && defined $yur)
    {
        for (my $j = 0; $j < $dx->numCols; $j++)
        {
            if (abs ($self->item(0, $j) - $xll) < $self->tol
            ||  abs ($self->item(0, $j) - $xur) < $self->tol)
            {
                $dx->set($j, 0);
            }
            
            if (abs ($self->item(1, $j) - $yll) < $self->tol
            ||  abs ($self->item(1, $j) - $yur) < $self->tol)
            {
                $dy->set($j, 0);
            }
        }
    }
    
    return $self->copy->plus($dx->vstack($dy));
}

################################################################################

# Min and max stuff.

sub minMax
{
    my $self = shift;
    
    my $min = $self->item(0, 0);
    my $iMin = 0;
    my $jMin = 0;
    
    my $max = $self->item(0, 0);
    my $iMax = 0;
    my $jMax = 0;
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            if ($self->item($i, $j) < $min)
            {
                $min = $self->item($i, $j);
                $iMin = $i;
                $jMin = $j;
            }
            
            if ($self->item($i, $j) > $max)
            {
                $max = $self->item($i, $j);
                $iMax = $i;
                $jMax = $j;
            }
        }
    }
    
    $$self{'min'}  = $min;
    $$self{'iMin'} = $iMin;
    $$self{'jMin'} = $jMin;
    
    $$self{'max'}  = $max;
    $$self{'iMax'} = $iMax;
    $$self{'jMax'} = $jMax;
}



sub min
{
    my $self = shift;
    return $$self{'min'} if defined $$self{'min'};
    $self->minMax;
    return $$self{'min'} if defined $$self{'min'};
    print "\nFailed to get the minimum value of the matrix.\n";  die;
}

sub iMin
{
    my $self = shift;
    return $$self{'iMin'} if defined $$self{'iMin'};
    $self->minMax;
    return $$self{'iMin'} if defined $$self{'iMin'};
    print "\nFailed to get the row index of the minimum value of the matrix.\n";  die;
}

sub jMin
{
    my $self = shift;
    return $$self{'jMin'} if defined $$self{'jMin'};
    $self->minMax;
    return $$self{'jMin'} if defined $$self{'jMin'};
    print "\nFailed to get the column index of the minimum value of the matrix.\n";  die;
}



sub max
{
    my $self = shift;
    return $$self{'max'} if defined $$self{'max'};
    $self->minMax;
    return $$self{'max'} if defined $$self{'max'};
    print "\nFailed to get the maximum value of the matrix.\n";  die;
}

sub iMax
{
    my $self = shift;
    return $$self{'iMax'} if defined $$self{'iMax'};
    $self->minMax;
    return $$self{'iMax'} if defined $$self{'iMax'};
    print "\nFailed to get the row index of the maximum value of the matrix.\n";  die;
}

sub jMax
{
    my $self = shift;
    return $$self{'jMax'} if defined $$self{'jMax'};
    $self->minMax;
    return $$self{'jMax'} if defined $$self{'jMax'};
    print "\nFailed to get the column index of the maximum value of the matrix.\n";  die;
}

################################################################################

sub toFile
{
    my $self = shift;
    my $filePath = shift;
    
    (open TXT, ">$filePath") || die;
    
    for (my $i = 0; $i < $self->numRows; $i++)
    {
        for (my $j = 0; $j < $self->numCols; $j++)
        {
            my $tmp = $self->item($i, $j);
            $tmp =~ s/\s+//g;
            print TXT $tmp;
            print TXT ' ' if $j != $self->numCols - 1;
        }
        
        print TXT "\n" if $i != $self->numRows - 1;
    }
    
    close TXT;
    sleep 1 while ! -e $filePath;
}

################################################################################

sub swapRows
{
    # Swap the rows of a matrix.
    my $self = shift;
    my $i = shift;
    my $j = shift;
    
    my $tmp = @{$self->items}[$i];
    @{$self->items}[$i] = @{$self->items}[$j];
    @{$self->items}[$j] = $tmp;
}



sub solve
{
    # Solve a linear system Ax = b for x.
    my $self = shift;
    my $rhs = shift;
    
    my $A = $self->copy;
    my $b = $rhs->copy;

    if ($A->numRows != $A->numCols || $A->numRows != $b->numRows) {
        print "Rows(A) = " . $A->numRows . "\n";
        print "Rows(b) = " . $b->numRows . "\n";
        print "\nA should be square.  Rows(A) should equal Rows(b).\n";  die;
    }
    
    my $nRows = $A->numRows;
    my $nCols = $A->numCols;

    # Apply row operations to transform A to upper triangular form.
    for (my $j = 0; $j < $nCols - 1; $j++)
    {
        # Find largest element in leftmost column, and make this the pivot row.
        # Apparently, this is important, because when I was not doing this, I
        # was getting noticeably different answers from python when the number
        # of subdomains was small.
        my $indMax = $j;
        for (my $i = $j+1; $i < $nRows; $i++)
        {
            if ((abs $A->item($i, $j)) > (abs $A->item($indMax, $j)))
            {
                $indMax = $i;
            }
        }
        # Swap row $indMax and row $j.
        $A->swapRows($indMax, $j);
        $b->swapRows($indMax, $j);
        # Zero out the $jth column.
        for (my $i = $j + 1; $i < $nRows; $i++)
        {
            my $factor = -1 * $A->item($i, $j) / $A->item($j, $j);
            for (my $k = $j; $k < $nCols; $k++)
            {
                $A->set($i, $k, $A->item($i, $k) + $factor * $A->item($j, $k));
            }
            for (my $ell = 0; $ell < $b->numCols; $ell++)
            {
                $b->set($i, $ell, $b->item($i, $ell) + $factor * $b->item($j, $ell));
            }
        }
    }

    # Use back substitution to finish solving for @x.
    my $x = $b;
    for (my $ell = 0; $ell < $b->numCols; $ell++)
    {
        $x->set($nRows-1, $ell, $b->item($nRows-1, $ell) / $A->item($nRows-1, $nRows-1));
    }
    my $k = $nRows;
    for (my $i = $nRows - 2; $i >= 0; $i--)
    {
        $k--;
        for (my $ell = 0; $ell < $b->numCols; $ell++)
        {
            my $dot = 0;
            for (my $j = $k; $j < $nRows; $j++)
            {
                $dot += ($A->item($i, $j) * $x->item($j, $ell));
            }
            $x->set($i, $ell, ($b->item($i, $ell) - $dot) / $A->item($i, $i));
        }
    }
    
    return $x;
}



sub leastSquares
{
    # "Solve" for x in an over-determined linear system $self * x = b.
    my $self = shift;
    my $b = shift;
    
    my $aT = $self->transpose;
    
    return $aT->dot($self)->solve($aT->dot($b));
}

################################################################################

return 1;
