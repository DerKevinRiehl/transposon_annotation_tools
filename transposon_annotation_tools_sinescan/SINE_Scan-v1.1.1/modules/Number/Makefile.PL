use 5.008;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile
    ( NAME         => 'Number::Format',
      VERSION_FROM => 'Format.pm', # finds $VERSION
      PREREQ_PM    => { 'Test::More' => 0,
                        'Carp'       => 0,
                        'POSIX'      => 0 },

      ($ExtUtils::MakeMaker::VERSION ge '6.31'?  (LICENSE => 'perl', ) : ()),

      AUTHOR       => 'William R. Ward <wrw@'.'cpan.org>',
      ABSTRACT     => 'Perl extension for formatting numbers',
    );
