This directory contains the source code for the TAL Effector Site Finder.

# License

All source code is available under an ISC license.

Copyright (c) 2011-2015, Daniel S. Standage <daniel.standage@gmail.com> and
Erin Doyle <edoyle@iastate.edu> with modifications by Nick Booher <njbooher@gmail.com>.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

# Description

TALESF is a C library for identifying potential binding sites for transcription
activator-like (TAL) effectors in a given genomic sequence. Compiling and
running this program requires a C compiler with OpenMP support (such as GCC 4.2
or higher). [libbcutils](https://github.com/boglab/cutils) must also be installed. To compile TALESF, enter the following commands in the directory containing this file.
```
  make
  make install
  make frontend
```

To install the Cython wrapper (only required if you want to use the Python frontend from [boglab_tools](https://github.com/boglab/talent_tools)), install [Cython](http://pypi.python.org/pypi/Cython) then run:

```
  cd cython_wrapper
  python setup.py build_ext
  python setup.py install
```

After compiling, you can run the program like this:
```
  ./pairedtalesf -o output_file genome_seq.fasta "HD NI NG NG NI HD NG NN NG NI NI NI NI N* NS N*" "N* NS N* NI NI NI NI NG NN NG HD NI NG NG NI HD"
```
This command will search the genome sequence saved in the file genome_seq.fasta for TALEN pair "HD NI NG NG NI HD NG NN NG NI NI NI NI N* NS N*" "N* NS N* NI NI NI NI NG NN NG HD NI NG NG NI HD" and save the output in files called output_file.txt and output_file.gff3. Note that the file genome_seq.fasta must be in [FASTA format](http://en.wikipedia.org/wiki/FASTA_format#Format).

For information on available options:
```
  ./pairedtalesf -h
```

By default, the output data is sorted by score. The 'sortfilter' script
is provided as a tool for (surprise!) sorting and filtering output from TALESF.
For a descriptive usage statement, enter the following command.
```
  ./sortfilter -h
```
