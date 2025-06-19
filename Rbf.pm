#!/mu/bin/perl

# A class to define an Rbf object, which is useful for working with RBFs.
# Greg Barnett
# December 2024

package Rbf;

use strict;
use warnings;

eval { use FindBin; };
use lib "$FindBin::Bin";
eval { use Matrix;  };

################################################################################

sub new
{
    my %self = ();
    my $self = \%self;
    $$self{'type'} = 'Rbf';
    bless $self;
    return $self;
}

################################################################################

sub test2d
{
    # /mu/bin/perl -e 'use lib "/home/gregorybarne/tmpCode"; use Rbf; Rbf::test2d();'
    # perl -e "use lib \"Z:\tmpCode\"; use Rbf; Rbf::test2d();"
    
    # Make a new Rbf object.
    my $phs = Rbf::new();

    # Set the number of dimensions (number of independent variables).
    $phs->setDimensions(2);

    # Set the exponent for the RBFs.  phi(r) = r**3 is exponent = 3, for example.
    $phs->setExponent(3);

    # Set the maximum polynomial degree that will be used in the approximation.
    $phs->setPolyDegree(1);

    # Bounding box for domain
    my $a = 0;
    my $b = 1;
    my $c = 0;
    my $d = 1;

    # Nodes
    my @x = ($a, $b, $b, $a, ($a + $b) / 2);
    my @y = ($c, $c, $d, $d, ($c + $d) / 2);
    $phs->setNodes(Matrix::new([\@x, \@y]));

    # Eval pts
    my $n = 5;                                    # number of columns (x-values)
    my $m = 5;                                    # number of rows (y-values)
    my $x = Matrix::linspace($a, $b, $n);
    my $y = Matrix::linspace($c, $d, $m);
    my ($xx, $yy) = Matrix::meshgrid($x, $y);
    $phs->setEvalPts($xx->flatten->vstack($yy->flatten));

    $phs->setFunctionValues(Rbf::firstTestFunction2d($phs->nodes));
    my $estimates = $phs->evaluate->reshape($m, $n);
    my $exact = Rbf::firstTestFunction2d($phs->evalPts)->reshape($m, $n);
    $estimates->round(8)->disp;
    $exact->round(8)->disp;
    $estimates->minus($exact)->round(8)->disp;
    
    print "--------------------------------------------\n\n";
    
    $phs->setFunctionValues(Rbf::secondTestFunction2d($phs->nodes));
    $estimates = $phs->evaluate->reshape($m, $n);
    $exact = Rbf::secondTestFunction2d($phs->evalPts)->reshape($m, $n);
    $estimates->round(8)->disp;
    $exact->round(8)->disp;
    $estimates->minus($exact)->round(8)->disp;
}

################################################################################

sub setDimensions
{
    if ((scalar @_) != 2)
    {
        print STDERR "This function requires exactly one input, the number of dimensions.";  die;
    }

    my $self = shift;
    my $dimensions = shift;

    if ($dimensions =~ /2|3|4|5|6|7|8|9/)
    {
        $$self{'dimensions'} = $dimensions;
    }
    else
    {
        print STDERR "Please choose 2, 3, 4, 5, 6, 7, 8, or 9 for the number of dimensions.";  die;
    }
}



sub dimensions
{
    my $self = shift;
    return $$self{'dimensions'} if defined $$self{'dimensions'};
    print STDERR "The number of dimensions has not yet been defined.";  die;
}

################################################################################

# Polynomial Methods

sub setPolyDegree
{
    if ((scalar @_) != 2)
    {
        print STDERR "This function requires exactly one input, the polynomial degree.";  die;
    }

    my $self = shift;
    my $polyDegree = shift;

    if ($polyDegree =~ /0|1/)
    {
        $$self{'polyDegree'} = $polyDegree;
    }
    else
    {
        print STDERR "Please choose either 0 or 1 for the polynomial degree.";  die;
    }
}



sub polyDegree
{
    my $self = shift;
    return $$self{'polyDegree'} if defined $$self{'polyDegree'};
    print STDERR "The polynomial degree has not yet been defined.";  die;
}



sub numPolys
{
    my $self = shift;
    return $$self{'numPolys'} if defined $$self{'numPolys'};
    $$self{'numPolys'} = 1                     if $self->polyDegree == 0;
    $$self{'numPolys'} = 1 + $self->dimensions if $self->polyDegree == 1;
    return $$self{'numPolys'};
}

################################################################################

# RBF Methods

sub setExponent
{
    if ((scalar @_) != 2)
    {
        print STDERR "This function requires exactly one input, the RBF exponent.";  die;
    }

    my $self = shift;
    my $exponent = shift;

    if ($exponent =~ /1|3/)
    {
        $$self{'exponent'} = $exponent;
    }
    else
    {
        print STDERR "Please choose either 1 or 3 for the RBF exponent.";  die;
    }
}



sub exponent
{
    my $self = shift;
    return $$self{'exponent'} if defined $$self{'exponent'};
    print STDERR "The RBF exponent has not yet been defined.";  die;
}



sub setNodes
{
    if ((scalar @_) != 2)
    {
        print STDERR "This method requires exactly one input, the nodes.  You gave " . (scalar @_ - 1);  die;
    }

    my $self = shift;
    my $nodes = shift;

    if ($$nodes{'type'} eq 'Matrix')
    {
        $$self{'nodes'} = $nodes;
    }
    elsif (ref $nodes)
    {
        $$self{'nodes'} = Matrix::new($nodes);
    }
    elsif (-e $nodes)
    {
        unless (open NOD, "<$nodes")
        {
            print STDERR "Failed to open nodes file \"$nodes\".";  die;
        }
        
        chomp (my @nod = <NOD>);
        close NOD;
        my @nodes = ();

        foreach my $line (@nod)
        {
            if ($line =~ /\S+/)
            {
                my @line = split /\s*\,\s*/, $line;
                @line = split /\s+/, $line if (scalar @line) == 1;
                push @nodes, \@line;
            }
        }

        $$self{'nodes'} = Matrix::new(\@nodes)->transpose;
    }
    else
    {
        print STDERR "Could not interpret input as either a matrix or a file path to nodes.";  die;
    }
}



sub nodes
{
    my $self = shift;
    return $$self{'nodes'} if defined $$self{'nodes'};
    print STDERR "The nodes have not yet been defined.";  die;
}

################################################################################

# Methods related to the evaluation points.

sub setEvalPts
{
    if ((scalar @_) != 2)
    {
        print STDERR "This method requires exactly one input, the evaluation points.  You gave " . (scalar @_ - 1);  die;
    }

    my $self = shift;
    my $evalPts = shift;

    if ($$evalPts{'type'} eq 'Matrix')
    {
        $$self{'evalPts'} = $evalPts;
    }
    elsif (ref $evalPts)
    {
        $$self{'evalPts'} = Matrix::new($evalPts);
    }
    elsif (-e $evalPts)
    {
        unless (open PTS, "<$evalPts")
        {
            print STDERR "Failed to open evaluation points file \"$evalPts\".";  die;
        }
        
        chomp (my @pts = <PTS>);
        close PTS;
        my @evalPts = ();

        foreach my $line (@pts)
        {
            if ($line =~ /\S+/)
            {
                my @line = split /\s*\,\s*/, $line;
                @line = split /\s+/, $line if (scalar @line) == 1;
                push @evalPts, \@line;
            }
        }

        $$self{'evalPts'} = Matrix::new(\@evalPts)->transpose;
    }
    else
    {
        print STDERR "Could not interpret input as either a matrix or a file path to evaluation points.";  die;
    }
}



sub evalPts
{
    my $self = shift;
    return $$self{'evalPts'} if defined $$self{'evalPts'};
    print STDERR "The evalPts have not yet been defined.";  die;
}

################################################################################

sub setFunctionValues
{
    if ((scalar @_) != 2)
    {
        print STDERR "This method requires exactly one input, the function values.  You gave " . (scalar @_ - 1);  die;
    }

    my $self = shift;
    my $functionValues = shift;

    if ($$functionValues{'type'} eq 'Matrix')
    {
        $$self{'functionValues'} = $functionValues;
    }
    elsif (ref $functionValues)
    {
        $$self{'functionValues'} = Matrix::new($functionValues);
    }
    elsif (-e $functionValues)
    {
        unless (open VAL, "<$functionValues")
        {
            print STDERR "Failed to open function values file \"$functionValues\".";  die;
        }
        
        chomp (my @functionValues = <VAL>);
        close VAL;
        my @f = ();

        foreach my $line (@functionValues)
        {
            if ($line =~ /\S+/)
            {
                push @f, $line;
            }
        }
        
        $$self{'functionValues'} = Matrix::new([\@f]);
    }
    else
    {
        print STDERR "Could not interpret input as either a matrix (row-vector) or a file path to function values.";  die;
    }
}



sub functionValues
{
    my $self = shift;
    return $$self{'functionValues'} if defined $$self{'functionValues'};
    print STDERR "The function values have not yet been defined.";  die;
}

################################################################################

sub phi
{
    # A simple scalar RBF function.
    if ((scalar @_) != 2)
    {
        print STDERR "phi function requires exactly one scalar input, a radius.  You gave " . (scalar @_ - 1);  die;
    }

    my $self = shift;
    my $r = shift;
    
    return $r ** $self->exponent;
}

################################################################################

sub polyMat
{
    # Function that returns the polynomial matrix at the input points.
    if ((scalar @_) != 2)
    {
        print STDERR "Function requires exactly one input, the points.  You gave " . (scalar @_ - 1);  die;
    }

    my $self = shift;
    my $pts = shift;

    if ($$pts{'type'} ne 'Matrix')
    {
        print STDERR "The type of the points input should be \"Matrix\".";  die;
    }

    my $p = Matrix::ones(1, $pts->numCols);

    if ($self->polyDegree >= 1)
    {
        for (my $k = 0; $k < $self->dimensions; $k++)
        {
            $p = $p->vstack($pts->rows([$k]));
        }
    }

    return $p;
}

################################################################################

sub aMat
{
    my $self = shift;
    return $$self{'aMat'} if defined $$self{'aMat'};

    my $n = $self->nodes->numCols;
    my $A = Matrix::alloc($n, $n);

    for (my $i = 0; $i < $n; $i++)
    {
        for (my $j = 0; $j < $n; $j++)
        {
            my $r2 = 0;

            for (my $k = 0; $k < $self->dimensions; $k++)
            {
                $r2 += ($self->nodes->item($k, $i) - $self->nodes->item($k, $j)) ** 2;
            }

            $A->set($i, $j, $self->phi(sqrt $r2));
        }
    }

    $A = $A->hstack($self->polyMat($self->nodes)->transpose);
    my $tmp = $self->polyMat($self->nodes)->hstack(Matrix::zeros($self->numPolys, $self->numPolys));
    $A = $A->vstack($tmp);

    $$self{'aMat'} = $A;
    return $A;
}

################################################################################

# sub coeffs
# {
    # my $self = shift;
    # return $$self{'coeffs'} if defined $$self{'coeffs'};

    # my $rhs = $self->functionValues->transpose->vstack(Matrix::zeros($self->numPolys, $self->functionValues->numRows));
    # $$self{'coeffs'} =  $self->aMat->solve($rhs);
    # return $$self{'coeffs'};
# }

################################################################################

sub bMat
{
    my $self = shift;
    return $$self{'bMat'} if defined $$self{'bMat'};

    my $B = Matrix::alloc($self->evalPts->numCols, $self->nodes->numCols);

    for (my $i = 0; $i < $self->evalPts->numCols; $i++)
    {
        for (my $j = 0; $j < $self->nodes->numCols; $j++)
        {
            my $r2 = 0;

            for (my $k = 0; $k < $self->dimensions; $k++)
            {
                $r2 += ($self->evalPts->item($k, $i) - $self->nodes->item($k, $j)) ** 2;
            }
            
            $B->set($i, $j, $self->phi(sqrt $r2));
        }
    }

    $B = $B->hstack($self->polyMat($self->evalPts)->transpose);

    $$self{'bMat'} = $B;
    return $B;
}

################################################################################

sub weights
{
    my $self = shift;
    return $$self{'w'} if defined $$self{'w'};
    
    $$self{'w'} = $self->aMat->solve($self->bMat->transpose)->transpose;
    return $$self{'w'};
}

################################################################################

sub evaluate
{
    my $self = shift;
    return $self->weights->dot($self->functionValues->transpose->vstack(Matrix::zeros($self->numPolys, 1)))
    # return $self->bMat->dot($self->coeffs)->transpose;
}

################################################################################

# STATIC SUBROUTINES

sub firstTestFunction2d
{
    if (scalar @_ != 1)
    {
        print STDERR "This function requires exactly one input, the nodes, given in matrix form.";  die;
    }

    my $pts = shift;

    return $pts->rows([0])->pow(2)->plus($pts->rows([1]));
}



sub secondTestFunction2d
{
    if (scalar @_ != 1)
    {
        print STDERR "This function requires exactly one input, the nodes, given in matrix form.";  die;
    }

    my $pts = shift;

    return $pts->rows([0])->cos->dotTimes($pts->rows([1])->sin);
}

################################################################################

return 1;


