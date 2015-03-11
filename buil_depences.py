perl_package_name = 'LWP::UserAgent'
toolshed_package_name = 'perl_lwp_useragent'
toolshed_package_version = '1.57'

!cpanm CPAN::FindDependencies
!export PERL5LIB=$HOME/perl5/lib/perl5

%%perl --out output --err error
use CPAN::FindDependencies;
@deps = CPAN::FindDependencies::finddeps('LWP::UserAgent');
foreach my $dep (@deps) {
    print ' ' x $dep->depth;
    print "http://www.cpan.org/authors/id/", $dep->distribution(), "\n";
}
from string import Template
template = Template("""
<tool_dependency>
    <package name="perl" version="5.18.1">
        <repository name="package_perl_5_18" owner="iuc" prior_installation_required="True" />
    </package>
    <package name="$toolshed_name" version="$toolshed_version">
        <install version="1.0">
            <actions>
                <action type="setup_perl_environment">
                    <repository name="package_perl_5_18" owner="iuc">
                        <package name="perl" version="5.18.1" />
                    </repository>
$package
                </action>
            </actions>
        </install>
        <readme><![CDATA[
            Perl package: $readme
        ]]>
        </readme>
    </package>
</tool_dependency>
""")
package_content = []
for line in output.split():
    line = line.strip()
    if line:
        package_content.append( '                    <package>%s</package>' % line )

content = dict(
            readme=perl_package_name, 
            package='\n'.join( reversed(package_content) ),
            toolshed_name=toolshed_package_name,
            toolshed_version=toolshed_package_version
        )
print template.substitute(content)
