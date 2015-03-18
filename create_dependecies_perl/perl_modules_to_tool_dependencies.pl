# usar esse link para rodar no ipython
# http://nbviewer.ipython.org/github/bgruening/notebooks/blob/master/Perl/perl_package_to_galaxy_tool_dependency_v2.ipynb

#!/usr/bin/env perl
use strict;
use warnings;
use CPAN::FindDependencies;

my $package_name = $ARGV[0];
if(! defined($package_name)){
    die "Please invoke with the command: \n\n\tperl $0 My::Module::Name;\n\n";
}

my @deps = CPAN::FindDependencies::finddeps($package_name, perl=>"5.18.1");
# Reverse ordering by depth, and mapped to the distribution url
my @ordered_deps = map {$_->distribution() } sort {$b->depth <=> $a->depth} @deps;

my $template = <<"EOL";
<tool_dependency>
    <package name="perl" version="5.18.1">
    <repository name="package_perl_5_18" owner="iuc" prior_installation_required="True" />
    </package>
    <package name="%s" version="%s">
    <install version="1.0">
    <actions>
    <action type="setup_perl_environment">
    <repository name="package_perl_5_18" owner="iuc">
    <package name="perl" version="5.18.1" />
    </repository>
    %s
    </action>
    </actions>
    </install>
    <readme><![CDATA[
		   Perl package: %s
	       ]]>
    </readme>
    </package>
    </tool_dependency>
    EOL

my $package_deps = join("\n", map{ " " x 20 . "<package>http://www.cpan.org/authors/id/" . $_ . "</package>"} @ordered_deps );

# Construct package name from passed value
my $package_pkgname = sprintf("perl_%s", lc($package_name));
$package_pkgname =~ s/::/_/g;

# splits P/PE/PEVANS/Scalar-List-Utils-1.41.tar.gz
my @tmp = split(/-/, $deps[0]->distribution);

# Grab 1.41.tar.gz
my $package_pkgversion = $tmp[-1];

# Strip .tar.gz ending
$package_pkgversion =~ s/.tar.gz//g;

printf($template, $package_pkgname, $package_pkgversion, $package_deps, $package_name);
