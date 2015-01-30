# Galaxy-Tools

Galaxy-Tools is a online repository for proteomics statistical analysis tools designed by LNBio's team on proteomics. The team is composed of:

  - Daniel Germano Travieso
  - Flavia Vischi Winck
  - Mateus Bellomo Ruivo

The team works on [CNPEM]'s [LNBio].

### Version
0.0.1

### Tech

Galaxy-Tools uses some sources to work:
* [Galaxy Project] - Genomics web-based framework for tools
* [R] - R Statistics analysis tool
* [Perl] - Perl Language
* [BioConductor] - Bioinformatics R Package
* Many More...

### Installation

You need mercurial and galaxy installed locally

```sh
$ hg clone https://bitbucket.org/galaxy/galaxy-dist/
```

```sh
$ cd galaxy-dist
$ hg update stable
$ sh run.sh
```

Put all the files under a *new-folder*
> /path/to/galaxy-dist/tools/*new-folder*

Change the file *tool_conf.xml.sample* to be:
'''xml
<?xml version='1.0' encoding='utf-8'?>
<toolbox>
  <section id="*new-folder*" name="Section Name">
    <tool file="*new-folder*/*tool-name*.xml" />
    [...]
  </section>
  [...]
</toolbox>
'''

### Todo's

 - Write Tests
 - Complete Statistical Tools
 - Generate Workflows

[CNPEM]:http://cnpem.br
[LNBio]:http://lnbio.cnpem.br/
[Galaxy Project]:http://galaxyproject.org/
[R]:http://www.r-project.org/
[Perl]:https://www.perl.org/
[BioConductor]:http://bioconductor.org/
