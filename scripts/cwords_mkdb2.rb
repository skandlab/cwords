##
## biological sequence background
## rank by length and GC-content ranks
##

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
options[:threads]=1
options[:bg_window]=100

$coptions = OptionParser.new do |opts|
  opts.banner = "Usage: cwords_mkdb [options]"
  opts.on("-w", "--wordsize ARG", "wordsize") { |o| options[:wordsize] = o.split(",").map{|x| x.to_i}.sort}
  opts.on("-s", "--seqfile ARG", "sequence file (FASTA format)") {|o| options[:seqfile] = o}
  opts.on("-t", "--threads ARG", "use multiple threads to parallelize computations") {|o| options[:threads] = o.to_i}
  opts.on("-i", "--window ARG", "sliding window size for control sequences") {|o| options[:bg_window] = o.to_i}
end

def show_help(msg="", code=0, io=STDOUT)
  io.puts "#{msg}\n#{$coptions}"
  exit(code)
end

$coptions.parse!(ARGV)
#mandatory parameters
[:seqfile].each{ |p| show_help("option '#{p}' mandatory") if options[p].nil?}

exit("seqfile must have fasta-format") if !options[:seqfile].match(/.fa$/)
bgw = options[:bg_window] #/2
dbdir = basedir + "/db/" + File.basename(options[:seqfile],'.fa') + "_bio#{options[:bg_window]}"
FileUtils.mkdir_p dbdir # create dir if it does not exist

###
### Main program
###

puts ">> Parameters"
options.each{|k,v| puts sprintf("%-20s:  %s",k,v) if !v.nil?}

# word id's
wids = Hash.new
begin
  wi = 0
  options[:wordsize].each do |ws|
    ['a','g','c','t'].rep_perm(ws) {|seqa| wids[seqa.join('')]=wi ; wi+=1 }
  end
end

# preprocess seqs, sort by length
seqs = Hash.new() {|h,k| h[k] = []}
IO.readlines(options[:seqfile],">")[1..-1].each do |s|
  ff = s.split("\n").map{|x| x.chomp}
  id = ff.shift
  seq = ff[0..-2].join('').downcase.gsub('u','t') # last field is ">"
  seq += ff[-1] if ff[-1] != '>' # check if last field is ">"
  # next if not nucleotide sequence, i.e. "unavailable"
  next if (seq.split('').uniq - ['a','c','g','t']).size > 0
  next if seq.size < 50 # lower bound
  seqs[seq] << id # remove duplicates
end

nseqs = seqs.size
puts "\n>> Computing statistics for #{nseqs} unique sequences"

puts "\n>> listing sequence neighbours"

# rank by length and GC-content ranks
seqs = seqs.to_a
length_ranks = seqs.map{|x| x[0].size}.ranks
comp_ranks = seqs.map{|x| x[0].split('').select{|c| 'gc'.include?(c)}.size.to_f/x[0].size}.ranks
sranks = (0..nseqs-1).to_a.map{|x| [length_ranks[x],comp_ranks[x],x]}

snbs = Array.new(nseqs) # for each seq, store ids of neighbours
pbar = ProgressBar.new("seqs",nseqs)
(0..nseqs-1).to_a.threach(options[:threads]) do |sidx|
  s = sranks[sidx]
  # first sorted neighbour is is s
  snbs[sidx] = sranks.sort_by{|x| (x[0] - s[0]).abs + (x[1] - s[1]).abs}[1,bgw+1].map{|x| x[2]}
  pbar.inc
end
pbar.finish

observed = Array.new(seqs.size) {Array.new(wids.size,0)}

puts "\n>> Enumerating words in sequences"
pbar = ProgressBar.new("seqs",nseqs)
(0..nseqs-1).to_a.threach(options[:threads]) do |sidx|
  seq = seqs[sidx][0]
  seqsize = seq.size
  options[:wordsize].each{|ws| (0..seqsize-ws).each{|pi| observed[sidx][wids[seq[pi, ws]]] += 1}}
  pbar.inc
end
pbar.finish

puts "\n>> Computing and storing statistics"
pbar = ProgressBar.new("words",wids.size)
wids.to_a.threach(options[:threads]) do |w,wid|
  File.open("#{dbdir}/#{w}.rnk", 'w') do |f|
    # for each seq, pick and compare to N closest seqs in
    # neighbourhood
    (0..nseqs-1).to_a.each do |sidx|
      nobs = observed[sidx][wid]
      next if nobs == 0
      ngte = snbs[sidx].select{|x| observed[x][wid] >= nobs}.size
      lods = -Math.log((1+ngte)/(bgw+1.0))
      seqs[sidx][1].each do |seqid|
        f.puts [seqid,nobs,lods].join("\t")
      end
    end
  end
  pbar.inc
end
pbar.finish
