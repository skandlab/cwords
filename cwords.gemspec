Gem::Specification.new do |s|
  s.name = 'cwords'
  s.version = '0.1.11'
  s.date = '2010-04-09'
  s.authors = ['Anders Jacobsen']
  s.email = 'andersmbj@gmail.com'
  s.summary = 'Word correlation analysis'
  s.description = 'Word correlation analysis in ranked nucleotide sequences (bioinformatics)'
  s.homepage = 'https://github.com/andersjacobsen/cWords/'
  s.files = ['README', 'LICENSE', 'bin/cwords', 'bin/cwords_mkdb', 'lib/ushuffle.jar', 'lib/wordRS-lib.rb', 'resources/genemap.tsv', 'resources/word_annotation.tsv', 'scripts/cwords.rb','scripts/cwords_mkdb.rb','scripts/cwords_mkdb_worker.rb','scripts/cluster_words.rb', 'scripts/complementary_words.rb']
  s.executables = ['cwords','cwords_mkdb']
  s.default_executable = ['bin/cwords']
  s.add_dependency('progressbar','>= 0.9.0')
end
