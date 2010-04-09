
srcdir = File.dirname(__FILE__)
basedir = srcdir + "/../"
libdir = basedir + 'lib/'
$LOAD_PATH << libdir

require 'wordRS-lib.rb'
require 'rubygems'
require 'progressbar'
require 'optparse'
require 'fileutils'
require 'java'
require libdir + 'ushuffle.jar'
java_import 'UShuffle'
us = UShuffle.new

###
### Main
###


#default options
options = Hash.new
options[:wordsize] = [7]
options[:seqfile] = nil
options[:partition] = nil
options[:stats] = ['p'] # p=p-value, z=z-score
options[:shuffles]=1000
options[:bg]=1 #mononucleotide shuffling

$coptions = OptionParser.new do |opts|
  opts.on("-w", "--wordsize ARG", "wordsize") { |o| options[:wordsize] = o.split(",").map{|x| x.to_i}.sort}
  opts.on("-s", "--seqfile ARG", "sequence file (FASTA format)") {|o| options[:seqfile] = o}
  opts.on("-p", "--partition ARG", "only process a partition (i.e. 5-10) of sequences") {|o| options[:partition] = o}
  opts.on("-a", "--stats ARG", "sequence file") {|o| options[:stats] = o.split(',')}
  opts.on("-u", "--shuffle ARG", "number of shuffles") {|o| options[:shuffles] = o.to_i}
  opts.on("-b", "--bg ARG", "background nucleotide model") {|o| options[:bg] = o.to_i}
end

def show_help(msg="", code=0, io=STDOUT)
  io.puts "#{msg}\n#{$coptions}"
  exit(code)
end

$coptions.parse!(ARGV)
#mandatory parameters
[:seqfile].each{ |p| show_help("option '#{p}' mandatory") if options[p].nil?}

exit("seqfile must have fasta-format") if !options[:seqfile].match(/.fa$/)
dbdir = basedir + "/db/" + File.basename(options[:seqfile],'.fa') + "_bg#{options[:bg]}"
FileUtils.mkdir_p dbdir # create dir if it does not exist

decimals = 6
bg = options[:bg]

# word id's
@wid = Hash.new
i = 0
options[:wordsize].each do |ws|
  ['a','g','c','t'].rep_perm(ws) {|seqa| @wid[seqa.join('')]=i ; i+=1 }
end

@seqs = IO.readlines(options[:seqfile],">")[1..-1]
if options[:partition]
  puts "partition #{options[:partition]}"
  (pstart,pstop) = options[:partition].split('-')
  @seqs = @seqs[pstart.to_i-1..pstop.to_i-1]
end

puts "computing statistics for #{@seqs.size} sequences"
pbar = ProgressBar.new("seqs",@seqs.size)

@seqs.each do |s|
  ff = s.split("\n")
  id = ff.shift
  seq = ff[0..-2].join('').downcase.gsub('u','t') # last field is ">"
  # next if not nucleotide sequence, i.e. "unavailable"
  next if (seq.split('').uniq - ['a','c','g','t']).size > 0
  
  #observed word counts
  @observed = Array.new(@wid.size,0)
  options[:wordsize].each{|ws| (0..seq.size-ws).each{|i| wid = @wid[seq[i, ws]]; @observed[wid] += 1 if not wid.nil?}}
  
  #expected word counts
  @expected = Array.new(@wid.size) {Array.new(options[:shuffles],0).to_statarray}
  us.init_shuffle(seq,bg)
  options[:shuffles].times do |si|
    seqsh = us.shuffle
    options[:wordsize].each{|ws| (0..seq.size-ws).each{|i| wid = @wid[seqsh[i, ws]]; @expected[wid][si] += 1 if not wid.nil?}}
  end
  
  #store results
  @wid.each do |w,wid|
   obs = @observed[wid]
   File.open("#{dbdir}/#{w}.rnk", 'a') {|f| f.puts [id,obs,@expected[wid].select{|x| x>=obs}.size,@expected[wid].to_statarray.mean].join("\t")}
  end
  pbar.inc
end
pbar.finish
