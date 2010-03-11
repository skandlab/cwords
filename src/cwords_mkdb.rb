#!/usr/bin/ruby

srcdir = File.dirname(__FILE__)
basedir = srcdir + "../"
libdir = basedir + 'lib/'
$LOAD_PATH << libdir

require 'wordRS-lib.rb'
require 'progressbar'
require 'optparse'
require 'fileutils'

tdir = basedir + '/tmp/'
FileUtils.mkdir_p tdir # create dir if it does not exist

###
### Main
###

#default options
options = Hash.new
options[:wordsize] = [7]
options[:seqfile] = nil
options[:partitions] = 1
options[:stats] = ['p'] # p=p
options[:ruby]='jruby --fast -J-Xmx1024m'
options[:shuffles]=5000
options[:bg]=1 #mononucleotide shuffling

$coptions = OptionParser.new do |opts|
  opts.on("-w", "--wordsize ARG", "wordsize") { |o| options[:wordsize] = o.split(",").map{|x| x.to_i}}
  opts.on("-s", "--seqfile ARG", "sequence file") {|o| options[:seqfile] = o}
  opts.on("-p", "--partitions ARG", "number of sequence partitions") {|o| options[:partitions] = o.to_i}
  opts.on("-a", "--stats ARG", "sequence file") {|o| options[:stats] = o.split('')}
  opts.on("-u", "--shuffle ARG", "number of shuffles") {|o| options[:shuffles] = o.to_i}
  opts.on("--ruby ARG", "ruby interpreter") {|o| options[:ruby] = o}
  opts.on("-b", "--bg ARG", "background nucleotide model") {|o| options[:bg] = o.to_i}
end

def show_help(msg="", code=0, io=STDOUT)
  io.puts "#{msg}\n#{$coptions}"
  exit(code)
end

$coptions.parse!(ARGV)
#mandatory parameters
[:seqfile].each{  |p| show_help("option '#{p}' mandatory") if options[p].nil?}

exit("seqfile must have fasta-format") if !options[:seqfile].match(/.fa$/)
dbname = File.basename(options[:seqfile],'.fa')
dbdir = basedir + "/db/" + dbname + "_bg#{options[:bg]}"
FileUtils.mkdir_p dbdir # create dir if it does not exist

n=options[:partitions]

# word id's
@seqs = IO.readlines(options[:seqfile],"\n>")
puts "#{@seqs.size} sequences"

puts "purging database ..."
options[:wordsize].each do |wordsize|
  ['a','g','c','t'].rep_perm(wordsize) {|seqa|  wf = "#{dbdir}/#{seqa.join('')}.rnk"; File.delete(wf) if File.exist?(wf)}
end
      
puts "starting #{n} processes ..."

cmd = "#{options[:ruby]} #{basedir}/scripts/wordsrus_mkdb.rb"
cmd += " -w #{options[:wordsize].join(',')} -s #{options[:seqfile]} -a #{options[:stats].join(",")} -u #{options[:shuffles]} --bg #{options[:bg]}"

stamp = Time.now.to_i

partsize = @seqs.size/n
cmds = []
(n-1).times do |i|
  cmds << cmd + " -p #{(i)*(partsize)+1}-#{(i+1)*(partsize)} &> #{tdir}#{dbname}_b#{options[:bg]}_#{i+1}_#{stamp}.dbout"
end
cmds << cmd + " -p #{partsize*(n-1)+1}-#{[n*(partsize),@seqs.size].max} &> #{tdir}#{dbname}_b#{options[:bg]}_#{n}_#{stamp}.dbout"
cmds.each do |c|
  p c
  exec c if fork.nil?
end

puts "Jobs started."
puts "Monitor with : tail #{tdir}#{dbname}_*#{stamp}.dbout"
