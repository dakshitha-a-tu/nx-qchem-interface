#!/usr/bin/env perl
#
#====================================================================================
#
# Interface between Q-Chem and NewtonX
#
# Adapted from the TPP interface template by Mario Barbatti, 2021.
#
# Authors:
# Dakshitha Abeygunewardane
# Spiridoula Matsika
#
# Matsika Lab,
# Temple University, USA.
# 2025.
#
#====================================================================================
#
use strict;
use warnings;
use File::Copy;
use File::Path qw(make_path);
use lib join( '/', $ENV{"NX"}, "lib" );
use colib_perl;

# We are using "strict." Thus, all variables must be declared.
my ( $mld, $mdle, $BASEDIR, $DEBUG, $JAD, $JND, $ctd, $fpar );
my ( $nt, $qchem_template, $qchem_inp, $qchem_out, $print_extracted );
my ( $nat, $istep, $nstat, $nstatdyn, $ndamp, $kt, $dt, $t, $tmax, $nintc );
my ( $nxrestart, $thres, $killstat, $timekill, $prog, $lvprt, $etot_jump, $etot_drift );
my ( @symb, @zn, @x, @y, @z, @Mass, $ia );
my ( $typeofinput, $type );
my ( @epot, $ns );
my ( $vdoth, $run_flag );
my ( $ncoup, $nc );
my ( @gx, @gy, @gz, @hx, @hy, @hz );
my ( $mem, $retval );
my (%progconf);

# Define basic variables
define_variables();

# Write debug messages
write_debug1();

# Read current geometry
read_geom();

# Prepare Q-Chem input
prepare_input();

# Run Q-Chem
run_program();

# Read and write energies
read_energy();

# Read and write gradients
read_gradients();

# Nonadiabatic couplings
treat_nacme();

#
#====================================================================================
#
#                               START SUBROUTINES
#
#====================================================================================
#

sub define_variables {
    #
    #====================================================================================
    #
    # Define basic variables to be used everywhere.
    #
    #------------------------------------------------------------------------------------
    #

    $nt              = 4; # number of threads for qchem
    $print_extracted = 1;  # Set to 1 to print extracted data, 0 to suppress
    $mld             = $ENV{"NX"};
    $mdle            = "run-qchem-eom.pl:";
    $qchem_template  = "JOB_NAD/qchem.inp";
    $qchem_inp       = "qchem.inp";
    $qchem_out       = "qchem.out";

    $BASEDIR = `pwd`;
    chomp($BASEDIR);

    $DEBUG = "DEBUG";
    $JAD   = "JOB_AD";
    $JND   = "JOB_NAD";
    $ctd   = "control.d";

    # Load current dynamics status (read control.d)
    ( $nat, $istep, $nstat, $nstatdyn, $ndamp, $kt, $dt, $t, $tmax, $nintc, $mem, $nxrestart, $thres, $killstat, $timekill, $prog, $lvprt, $etot_jump, $etot_drift ) = load_status( $ctd, $mdle );

    %progconf = prog_config($prog);
    $fpar     = $progconf{parfile};
}

sub write_debug1 {
    #
    #====================================================================================
    #
    # Write debug messages.
    #
    #------------------------------------------------------------------------------------
    #
    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle has taken over\n",      $istep, $kt ); }
    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle running in $BASEDIR\n", $istep, $kt ); }
    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle beginning here\n",      $istep, $kt ); }
    if ( $lvprt >= 3 ) {
        print_STDOUT( "$mdle Dynamics control: \n",                                       $istep, $kt );
        print_STDOUT( "$nat $istep $nstat $nstatdyn $ndamp $kt $dt $t $tmax \n",          $istep, $kt );
        print_STDOUT( "$nintc $mem $nxrestart $thres $killstat $timekill $prog $lvprt\n", $istep, $kt );
        print_STDOUT( "$etot_jump $etot_drift\n\n",                                       $istep, $kt );
    }
}

sub read_geom {
    #
    #====================================================================================
    # This subroutine reads the current geometry and populates the
    # following vectors:
    # @symb      - atomic symbol
    # @zn        - atomic number
    # @x, @y, @z - cartesian coordinates (au)
    # @Mass      - atomic mass (amu)
    # The vector index runs from $ia = 0 to $ia = $nat-1.
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Reading current geometry.\n", $istep, $kt ); }

    open( GM, "geom" ) or die "$mdle Cannot open geom";
    $ia = 0;
    while (<GM>) {
        chomp;
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        ( $symb[$ia], $zn[$ia], $x[$ia], $y[$ia], $z[$ia], $Mass[$ia] ) = split( /\s+/, $_ );
        $ia++;
    }
    close(GM);
}

sub prepare_input {
    #
    #====================================================================================
    # This subroutine prepares the Q-Chem input by modifying a template.
    # It updates the geometry and the state for which properties are calculated.
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Preparing input for Q-Chem.\n", $istep, $kt ); }

    # Check if a template file exists
    die "$mdle Q-Chem template file '$qchem_template' not found!" unless -e $qchem_template;

    open( my $fh_tpl, '<', $qchem_template ) or die "Could not open '$qchem_template': $!";
    my @template_lines = <$fh_tpl>;
    close $fh_tpl;

    open( my $fh_inp, '>', $qchem_inp ) or die "Could not create '$qchem_inp': $!";

    my $au2ang = units("au2ang");
    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Unit conversion factor (au2ang): $au2ang\n", $istep, $kt ); }
    my $in_molecule_section = 0; # Flag to know if we are inside the molecule block

    # Find the charge and multiplicity line from the template first
    my $charge_mult_line = "";
    for (my $i = 0; $i < @template_lines; $i++) {
        if ($template_lines[$i] =~ /\$molecule/i) {
            $charge_mult_line = $template_lines[$i+1]; # The next line should be charge/mult
            last;
        }
    }
    die "Could not find charge/multiplicity line in template" unless ($charge_mult_line ne "");


    foreach my $line (@template_lines) {
        if ( $line =~ /\$molecule/i ) {
            $in_molecule_section = 1;
            print $fh_inp $line; # Print the $molecule line
            print $fh_inp $charge_mult_line; # Print the charge/multiplicity line

            # Write the new geometry
            for my $i ( 0 .. $#symb ) {
                my $x_ang = $x[$i] * $au2ang;
                my $y_ang = $y[$i] * $au2ang;
                my $z_ang = $z[$i] * $au2ang;
                printf $fh_inp "%-4s %15.8f %15.8f %15.8f\n", $symb[$i], $x_ang, $y_ang, $z_ang;
            }
            next; # Move to the next line in the template
        }

        if ( $line =~ /\$end/i && $in_molecule_section ) {
            $in_molecule_section = 0;
            print $fh_inp $line; # Print the $end line
            next;
        }

        if ($in_molecule_section) {
            # This will skip the old charge/mult line and the old geometry lines
            next;
        }

        # Handle other lines outside the molecule block
        if ( $line =~ /CC_STATE_TO_OPT/i ) {
            print $fh_inp "CC_STATE_TO_OPT   [1,$nstatdyn]\n";
        }
        else {
            print $fh_inp $line;
        }
    }

    close $fh_inp;
}

sub find_typeofinput {
    #
    #====================================================================================
    #
    # Check type of dynamics.
    #
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Checking type of dynamics.\n", $istep, $kt ); }

    # Read $type from type_of_dyn.out
    my $typeout = "type_of_dyn.out";
    $type = 2;
    if ( -s $typeout ) {
        open( INP, $typeout ) or die "$mdle Cannot open $typeout";
        $_ = <INP>;
        chomp;
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        my $type = $_;
        close(INP);
    }

    # Define type of input
    if ( ( -e $JAD ) and ( !-e $JND ) ) {
        $typeofinput = "jad";
    }
    elsif ( ( !-e $JAD ) and ( -e $JND ) ) {
        $typeofinput = "jnd";
    }
    elsif ( ( -e $JAD ) and ( -e $JND ) ) {
        if ( $type == 1 ) {
            $typeofinput = "jad";
        }
        elsif ( $type != 1 ) {
            $typeofinput = "jnd";
        }
    }
}

sub run_program {
    #
    #====================================================================================
    #
    # This subroutine executes Q-Chem and archives the output.
    #
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Executing Q-Chem.\n", $istep, $kt ); }

    find_typeofinput();

    if ( $typeofinput eq "jad" ) {
        die "The interface does not handle adiabatic dynamics: $!";
    }

    # Run Q-Chem
    $retval = callprog("", "qchem -nt $nt $qchem_inp $qchem_out > qchem.status", $mdle );
    #$retval = callprog("", "qchem $qchem_inp $qchem_out > qchem.status", $mdle );

    if ( $retval != 0 ) {
        print_STDOUT( "Error in the Q-Chem execution!\n", $istep, $kt );
        die;
    }

    # Archive the output file from the current step
    my $archive_dir = "../INFO_RESTART/qchem_outputs";
    make_path($archive_dir) unless -d $archive_dir;
    my $dest_file = "$archive_dir/$t.out";
    copy($qchem_out, $dest_file) or die "Failed to copy $qchem_out to $dest_file: $!";
    if ($lvprt >= 3) {
        print_STDOUT("$mdle Archived current output to $dest_file\n", $istep, $kt);
    }
}

sub read_energy {
    #
    #====================================================================================
    #
    # This subroutine reads the Q-Chem output to get the potential energies.
    # Energies must be in Hartree.
    #
    #------------------------------------------------------------------------------------
    #
    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Reading Epot.\n", $istep, $kt ); }

    open( my $fh_out, '<', $qchem_out ) or die "Could not open '$qchem_out': $!";
    my @lines = <$fh_out>;
    close $fh_out;

    my $line_num = 0;
    my $state_idx = 0;
    foreach my $line (@lines) {
        $line_num++;
        if ( $line =~ /(?:EOMIP|EOMEE) transition\s+\d+\/A/ ) {
            my $energy_line = $lines[$line_num]; # The energy is on the next line
            if ( $energy_line =~ /Total energy\s+=\s+(-?\d+\.\d+)\s+a\.u\./ ) {
                $epot[$state_idx] = $1;
                $state_idx++;
            }
        }
    }

    if ( $state_idx != $nstat ) {
        die "$mdle Expected $nstat states but found $state_idx energies in $qchem_out";
    }

    write_energy();
}

sub write_energy {
    #
    #====================================================================================
    #
    # This subroutine writes the potential energies to epot file and prints them.
    #
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Writing Epot.\n", $istep, $kt ); }

    my $file = "epot";
    if ( -s $file ) {
        copy( $file, "oldepot" ) or die "Copy failed: $!";
    }
    open( OUT, ">$file" ) or die "Cannot write to $file";
    foreach (@epot) {
        print OUT "$_\n";
    }
    close(OUT);
    copy( $file, "newepot" ) or die "Copy failed: $!";

    # Print energies to STDOUT if requested
    if ($print_extracted) {
        print_STDOUT( "State energies (a.u.):\n", $istep, $kt );
        for (my $i=0; $i < @epot; $i++) {
            my $state_num = $i + 1;
            my $line_to_print = sprintf("  State %d: %s\n", $state_num, $epot[$i]);
            print_STDOUT( $line_to_print, $istep, $kt );
        }
        print_STDOUT("\n", $istep, $kt);
    }
}

sub read_gradients {
    #
    #====================================================================================
    #
    # This subroutine reads the energy gradients from the Q-Chem output.
    # It parses the "G_I, a.u." and "G_J, a.u." blocks associated with NAC calculations.
    #
    #------------------------------------------------------------------------------------
    #
    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Reading Gradients.\n", $istep, $kt ); }

    # Initialize array with zeros
    for ( $ns = 0; $ns <= $nstat - 1; $ns++ ) {
        for ( $ia = 0; $ia <= $nat - 1; $ia++ ) {
            $gx[$ns][$ia] = 0.0;
            $gy[$ns][$ia] = 0.0;
            $gz[$ns][$ia] = 0.0;
        }
    }

    open( my $fh_out, '<', $qchem_out ) or die "Could not open '$qchem_out': $!";
    my @lines = <$fh_out>;
    close $fh_out;

    my $curr_state_A = -1;
    my $curr_state_B = -1;
    my $line_num = 0;
    my $found_gradients = 0;

    while ($line_num < @lines) {
        my $line = $lines[$line_num];

        # Capture State A and State B indices
        # Regex matches lines like: "State A: eomip_ccsd/a: 3/A"
        if ( $line =~ /State A:.*:\s+(\d+)\// ) {
            $curr_state_A = $1;
        }
        if ( $line =~ /State B:.*:\s+(\d+)\// ) {
            $curr_state_B = $1;
        }

        # Check for State I Gradient Block
        if ( $line =~ /^\s*G_I,\s*a\.u\./ ) {
            if ($curr_state_A == -1) {
                die "$mdle Found G_I gradient block but State A is undefined!";
            }
            $line_num += 3; # Skip Header, Atom Label, and Separator lines

            for ( $ia = 0; $ia < $nat; $ia++ ) {
                my $dline = $lines[$line_num + $ia];
                # Split on whitespace.
                my @F = split( /\s+/, $dline );
                # Format: Atom_Index X Y Z. F[0] might be empty if leading space.
                my $idx_offset = ($F[0] eq "") ? 1 : 0;

                $gx[ $curr_state_A - 1 ][$ia] = $F[$idx_offset + 1];
                $gy[ $curr_state_A - 1 ][$ia] = $F[$idx_offset + 2];
                $gz[ $curr_state_A - 1 ][$ia] = $F[$idx_offset + 3];
            }
            $found_gradients = 1;
        }

        # Check for State J Gradient Block
        elsif ( $line =~ /^\s*G_J,\s*a\.u\./ ) {
            if ($curr_state_B == -1) {
                die "$mdle Found G_J gradient block but State B is undefined!";
            }
            $line_num += 3; # Skip Header, Atom Label, and Separator lines

            for ( $ia = 0; $ia < $nat; $ia++ ) {
                my $dline = $lines[$line_num + $ia];
                my @F = split( /\s+/, $dline );
                my $idx_offset = ($F[0] eq "") ? 1 : 0;

                $gx[ $curr_state_B - 1 ][$ia] = $F[$idx_offset + 1];
                $gy[ $curr_state_B - 1 ][$ia] = $F[$idx_offset + 2];
                $gz[ $curr_state_B - 1 ][$ia] = $F[$idx_offset + 3];
            }
            $found_gradients = 1;
        }

        $line_num++;
    }

    # Check if we are running a single state calculation (non-NAC) where only 'Final gradient' exists
    # If we haven't found any G_I/G_J blocks, look for a generic Final gradient for the target state.
    if (!$found_gradients) {
        if ( $lvprt >= 3 ) { print_STDOUT( "$mdle No G_I/G_J blocks found. Searching for 'Final gradient' (Single state assumption).\n", $istep, $kt ); }
        $line_num = 0;
        while ($line_num < @lines) {
            if ( $lines[$line_num] =~ /Final gradient/ ) {
                # In standard output, Final gradient usually refers to the optimized state.
                # We assume this is nstatdyn.
                $line_num += 3;
                for ( $ia = 0; $ia < $nat; $ia++ ) {
                    my $dline = $lines[$line_num + $ia];
                    my @F = split( /\s+/, $dline );
                    my $idx_offset = ($F[0] eq "") ? 1 : 0;
                    $gx[ $nstatdyn - 1 ][$ia] = $F[$idx_offset + 1];
                    $gy[ $nstatdyn - 1 ][$ia] = $F[$idx_offset + 2];
                    $gz[ $nstatdyn - 1 ][$ia] = $F[$idx_offset + 3];
                }
                $found_gradients = 1;
                last;
            }
            $line_num++;
        }
    }

    die "$mdle Could not find any gradient blocks in $qchem_out" unless $found_gradients;

    write_grad();
}

sub write_grad {
    #
    #====================================================================================
    #
    # This subroutine writes the energy gradients to grad and grad.all files and prints them.
    #
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Writing Gradients.\n", $istep, $kt ); }

    open( OUT1, ">grad.all" ) or die "Cannot write to grad.all";
    open( OUT2, ">grad" )     or die "Cannot write to grad";

    for ( $ns = 0; $ns <= $nstat - 1; $ns++ ) {
        for ( $ia = 0; $ia <= $nat - 1; $ia++ ) {
            print OUT1 "$gx[$ns][$ia]  $gy[$ns][$ia]  $gz[$ns][$ia]\n";
            if ( $ns == $nstatdyn - 1 ) {
                print OUT2 "$gx[$ns][$ia]  $gy[$ns][$ia]  $gz[$ns][$ia]\n";
            }
        }
    }

    close(OUT1);
    close(OUT2);

    # Print gradient of current state to STDOUT if requested
    if ($print_extracted) {
        print_STDOUT( "Gradient of state $nstatdyn (a.u.):\n", $istep, $kt );
        print_STDOUT( "Atom          X               Y               Z\n", $istep, $kt );
        print_STDOUT( "--------------------------------------------------\n", $istep, $kt );
        for ( $ia = 0; $ia <= $nat - 1; $ia++ ) {
            my $atom_num = $ia + 1;
            my $line_to_print = sprintf("%-4d %15.8f %15.8f %15.8f\n", $atom_num, $gx[$nstatdyn-1][$ia], $gy[$nstatdyn-1][$ia], $gz[$nstatdyn-1][$ia]);
            print_STDOUT($line_to_print, $istep, $kt);
        }
        print_STDOUT("\n", $istep, $kt);
    }
}

sub treat_nacme {
    #
    #====================================================================================
    #
    # This subroutine processes the NACME.
    #
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Treating couplings.\n", $istep, $kt ); }

    $run_flag = read_single_value( "which_run_is_that", "" );

    if ( $run_flag ne "second run" ) {
        $ncoup = $nstat * ( $nstat - 1 ) / 2;

        # Intialize arrays
        for ( $nc = 0; $nc <= $ncoup - 1; $nc++ ) {
            for ( $ia = 0; $ia <= $nat - 1; $ia++ ) {
                $hx[$nc][$ia] = 0.0;
                $hy[$nc][$ia] = 0.0;
                $hz[$nc][$ia] = 0.0;
            }
        }

        if ( -s "nad_vectors" ) {
            copy( "nad_vectors", "oldh" ) or die "Copy of nad_vectors to oldh failed: $!";
        }
        elsif ( !-s "nad_vectors" ) {
            write_nacme("oldh");
        }

        read_nacme();

        write_nacme("nad_vectors");
        adjust_phase();
        copy( "nad_vectors", "newh" ) or die "Copy failed: $!";
    }
}

sub read_nacme {
    #
    #====================================================================================
    #
    # This subroutine reads all available nonadiabatic coupling vectors from Q-Chem output.
    # It parses "State A" / "State B" sections and the "NAC d^x_IJ (CI part)" blocks.
    #
    #------------------------------------------------------------------------------------
    #
    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Reading coupling vectors.\n", $istep, $kt ); }

    open( my $fh_out, '<', $qchem_out ) or die "Could not open '$qchem_out': $!";
    my @lines = <$fh_out>;
    close $fh_out;

    my $line_num = 0;
    my $curr_state_A = -1;
    my $curr_state_B = -1;

    while ($line_num < @lines) {
        my $line = $lines[$line_num];

        # Capture State A and State B indices
        if ( $line =~ /State A:.*:\s+(\d+)\// ) {
            $curr_state_A = $1;
        }
        if ( $line =~ /State B:.*:\s+(\d+)\// ) {
            $curr_state_B = $1;
        }

        # Find the NAC table within this block
        if ( $line =~ /^\s*NAC d\^x_IJ \(CI part\), a\.u\./ ) {

            if ($curr_state_A == -1 || $curr_state_B == -1) {
                 if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Found NAC block but states undefined. Skipping.\n", $istep, $kt ); }
                 $line_num++;
                 next;
            }

            $line_num += 3; # Skip Header, Atom Label, and Separator lines

            my $higher_state = ($curr_state_A > $curr_state_B) ? $curr_state_A : $curr_state_B;
            my $lower_state = ($curr_state_A < $curr_state_B) ? $curr_state_A : $curr_state_B;

            # NewtonX uses lower triangular matrix indexing
            # nc = (i-1)*i/2 + j - 1 for i > j
            my $nc_index = (($higher_state-1)*($higher_state-2)/2) + $lower_state -1;

            for (my $ia = 0; $ia < $nat; $ia++) {
                my $dline = $lines[$line_num + $ia];
                my @F = split(/\s+/, $dline);
                my $idx_offset = ($F[0] eq "") ? 1 : 0;

                # QChem computes d_IJ = <I|grad|J>
                # NewtonX wants h_IJ = <I|grad|J> (usually) but checks phases carefully.
                # Standard convention for NX: h_IJ = <I|d/dR|J>
                # We store h_IJ for I > J.

                # Note: d_JI = -d_IJ.
                # If Q-Chem printed d_AB where A > B, we use it directly.
                # If Q-Chem printed d_AB where A < B, then d_BA = -d_AB.

                # Check Q-Chem print order vs State ID.
                # The output says "NAC d^x_IJ". Usually I=StateA, J=StateB.
                # If StateA=3, StateB=1 (A>B), we read d_31.
                # If StateA=1, StateB=3 (A<B), we read d_13.

                # We need h_31.
                # h_31 = d_31.
                # h_31 = -d_13.

                my $val_x = $F[$idx_offset + 1];
                my $val_y = $F[$idx_offset + 2];
                my $val_z = $F[$idx_offset + 3];

                my $sign = 1.0;
                if ($curr_state_A < $curr_state_B) {
                    $sign = -1.0;
                }

                $hx[$nc_index][$ia] = $sign * $val_x;
                $hy[$nc_index][$ia] = $sign * $val_y;
                $hz[$nc_index][$ia] = $sign * $val_z;
            }
        }
        $line_num++;
    }
}

sub write_nacme {
    #
    #====================================================================================
    #
    # This subroutine writes @h to a file and prints the non-zero couplings.
    #
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Writing coupling vectors.\n", $istep, $kt ); }

    my ($file) = @_;
    open( OUT, ">$file" ) or die "Cannot write to $file";

    for ( $nc = 0; $nc <= $ncoup - 1; $nc++ ) {
        for ( $ia = 0; $ia <= $nat - 1; $ia++ ) {
            print OUT "$hx[$nc][$ia]  $hy[$nc][$ia]  $hz[$nc][$ia]\n";
        }
    }
    close(OUT);

    # Print non-zero NACs to STDOUT if requested
    # Only print for the final "nad_vectors" file to avoid duplication during initialization
    if ($print_extracted && $file eq "nad_vectors") {
        print_STDOUT( "Non-adiabatic couplings (a.u.):\n", $istep, $kt );
        my $printed_header = 0;
        my $nc_index = 0;
        for (my $i = 2; $i <= $nstat; $i++) {
            for (my $j = 1; $j < $i; $j++) {

                # Filter: Only process if this pair involves the current dynamics state (nstatdyn)
                # But we still need to increment nc_index for every pair!
                my $is_relevant = ($i == $nstatdyn || $j == $nstatdyn);

                if ($is_relevant) {
                    my $is_zero = 1;
                    for ($ia = 0; $ia <= $nat - 1; $ia++) {
                        if ($hx[$nc_index][$ia] != 0.0 || $hy[$nc_index][$ia] != 0.0 || $hz[$nc_index][$ia] != 0.0) {
                            $is_zero = 0;
                            last;
                        }
                    }

                    if (!$is_zero) {
                        if (!$printed_header) {
                            print_STDOUT( "--------------------------------------------------\n", $istep, $kt );
                            $printed_header = 1;
                        }
                        my $header = sprintf("Coupling between state %d and state %d:\n", $j, $i);
                        print_STDOUT( $header, $istep, $kt );
                        print_STDOUT( "Atom          X               Y               Z\n", $istep, $kt );
                        for ($ia = 0; $ia <= $nat - 1; $ia++) {
                            my $atom_num = $ia + 1;
                            my $line_to_print = sprintf("%-4d %15.8f %15.8f %15.8f\n", $atom_num, $hx[$nc_index][$ia], $hy[$nc_index][$ia], $hz[$nc_index][$ia]);
                            print_STDOUT($line_to_print, $istep, $kt);
                        }
                        print_STDOUT( "--------------------------------------------------\n", $istep, $kt );
                    }
                }
                $nc_index++;
            }
        }
        print_STDOUT("\n", $istep, $kt);
    }
}

sub adjust_phase {
    #
    #====================================================================================
    #
    # This subroutine computes cos(q) where q is the angle between couplings at t and t-Dt.
    #
    #------------------------------------------------------------------------------------
    #

    if ( $lvprt >= 3 ) { print_STDOUT( "$mdle Fixing the coupling phase.\n", $istep, $kt ); }

    copy( "nad_vectors", "newh" ) or die "Copy failed: $!";

    $retval = callprogsp( $mld, "escalar", $mdle );
    if ( $retval != 0 ) {
        die "$mdle is dying now\n";
    }
    copy( "nadv", "nad_vectors" ) or die "Copy failed: $!";

    if ( $lvprt >= 3 ) { copy( "escalar.log", "../$DEBUG/." ) or die "Copy failed: $!"; }
}