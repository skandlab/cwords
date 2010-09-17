
# dwords - Discriminatory Words
# cwords - Correlating Words

#!/usr/bin/env jruby

srcdir = File.dirname(__FILE__)
basedir = srcdir + "/../"
libdir = basedir + '/lib/'
$LOAD_PATH << libdir

require 'wordRS-lib.rb'
require 'rubygems'
require 'progressbar'
require 'optparse'
require 'java'
require libdir + 'ushuffle.jar'
java_import 'UShuffle'

#default options
options = Hash.new
options[:wordsize] = [7]
options[:split_words]=nil
options[:scoring_scheme] = 'obs'
options[:bg_samples]=1000
options[:bg]=0 # 0=bio, 1,2,3, ...
options[:geneset]=nil
options[:seqfile]=nil
options[:report_words]=nil
options[:plot_words]=nil
options[:dump]=nil
options[:threads]=1

$coptions = OptionParser.new do |opts|
  opts.banner = "Usage: cwords [options]"
  
  # analysis settings
  opts.on("-c", "--scoring_scheme ARG", "scoring scheme") {|o| options[:scoring_scheme] = o} 
  opts.on("-p", "--bg_samples ARG", "samples in background distribution") {|o| options[:bg_samples] = o.to_i}
  opts.on("-b", "--bg ARG", "background model") {|o| options[:bg] = o.to_i}

  opts.on("-w", "--wordsize ARG", "wordsize") { |o| options[:wordsize] = o.split(",").map{|x| x.to_i}}
  opts.on("-t", "--threads ARG", "use multiple threads to parallelize computations") {|o| options[:threads] = o.to_i}
  
  # files and directories
  opts.on("-g", "--geneset ARG", "gene set") {|o| options[:geneset] = o}
  opts.on("-s", "--seqfile ARG", "sequence file") {|o| options[:seqfile] = o}
  opts.on("-d", "--db ARG", "word database") {|o| options[:db] = o}
  
  # output control
  opts.on("-u", "--dump ARG", "dump top words") {|o| options[:dump] = o.to_i}
  opts.on(      "--report_words ARG", "report on words (comma separated)") {|o| options[:report_words] = o.split(',')}
end

def show_help(msg="", code=0, io=STDOUT)
  io.puts "#{msg}\n#{$coptions}"
  exit(code)
end

begin
  $coptions.parse!(ARGV)
rescue OptionParser::ParseError => error
  puts error.message
  puts $coptions
  exit
end

# mandatory parameters
[:seqfile].each{|p| show_help("option '#{p}' mandatory") if options[p].nil?}
show_help("scoring scheme must be one of: obs,bin") if !(['obs','bin'].include?(options[:scoring_scheme]))

# if bg = 0, require geneset
show_help("geneset file required for bg=0") if options[:bg] == 0 and !options[:geneset]

# get filename without directory
rankfilename = File.basename(options[:geneset] || options[:seqfile])

output_top = 10
annofile = basedir + "/resources/" + "word_annotation.tsv" #annotation
tidfile = basedir + "/resources/" + "genemap.tsv"
sequences = nil
nwords = options[:wordsize].map{|x| 4**x}.to_statarray.sum
threads=options[:threads]
bg_samples = options[:bg_samples]

###
### Main program
###

puts ">> Parameters"
options.each{|k,v| puts sprintf("%-20s:  %s",k,v) if !v.nil?}

# read in mirbase seed family
word_annotation = Hash.new("") # seq => family
IO.readlines(annofile).each{|l| word_annotation[l.split("\t")[0]] = l.split("\t")[1]}

# initialize word id hash, word sequence => word id (0..nwords-1)
wids = Hash.new
begin
  wi = 0
  options[:wordsize].each{|ws| ['a','g','c','t'].rep_perm(ws) {|seqa| wids[seqa.join('')]=wi ; wi+=1 }}
end

puts "\n>> Reading sequences ..."
sequences = Hash.new # id => seq
IO.readlines(options[:seqfile],">")[1..-1].each do  |s|
  ff = s.split("\n").map{|x| x.chomp}
  id = ff.shift
  seq = ff[0..-2].join('').downcase.gsub('u','t')
  seq += ff[-1] if ff[-1] != '>' # check if last field is ">"
  # next if not nucleotide sequence, i.e. "unavailable"
  next if (seq.split('').uniq - ['a','c','g','t']).size > 0
  # hash ensures sequence ids are unique
  sequences[id] = seq
end  
nseqs = sequences.size
seqs = sequences.to_a # id , seq

zscores = Array.new(nwords)

if options[:bg] == 0 # bio bg
  puts "\n>> Using similar sequences as background ..."
  
  geneset = []
  begin
    # build internal id hash map
    external_ids = Hash.new
    IO.readlines(tidfile).each do |l|
      tid = l.split(" ")[0]
      l.chomp.split(" ")[1].split(",").each{|extid| external_ids[extid] = tid}
    end
    internal_ids = external_ids.invert # allowed internal ids
    
    # sequence ids
    seqids = Hash.new # id => idx
    seqs.each_with_index{|x,idx| seqids[x[0]] = idx}
    
    # geneset ids, try mapping directly to sequence ids
    # then try extid map hash, compact removes non-mappable entries
    genesetf = IO.readlines(options[:geneset]).map{|l| l.split("\t").first.chomp}
    geneset = genesetf.map{|x| seqids[x] || seqids[external_ids[x]]}.compact
    
    puts "mapped #{geneset.size} (of #{genesetf.size}) IDs in geneset to sequences"
  end
  
  # rank by length and GC-content ranks
  puts "\n>> Listing sequence neighbours ..."
  seqnbs = Array.new(geneset.size) # for each geneset seq, store ids of neighbours
  begin
    length_ranks = seqs.map{|x| x[1].size}.ranks
    comp_ranks = seqs.map{|x| x[1].split('').select{|c| 'gc'.include?(c)}.size.to_f/x[1].size}.ranks
    sranks = (0..nseqs-1).to_a.map{|idx| [length_ranks[idx],comp_ranks[idx],idx]}
    # remove geneset from background
    bgranks = sranks.select{|x| !geneset.include?(x[2])}
    pbar = ProgressBar.new("progress",geneset.size)
    (0..geneset.size-1).to_a.threach(threads) do |sidx|
      s = sranks[sidx]
      nbs = bgranks.sort_by{|x| (x[0] - s[0]).abs + (x[1] - s[1]).abs}[0,200]
      # shuffle to randomly group neighbours
      seqnbs[sidx] = (1..bg_samples).to_a.map{|bgi| nbs[rand(nbs.size-1)][2]}
      pbar.inc
    end
    pbar.finish
  end

  puts "\n>> Enumerating words in sequences"
  wordobs = Array.new(nseqs) {Array.new(nwords,0)}
  begin
  pbar = ProgressBar.new("progress",nseqs)
    (0..nseqs-1).to_a.each do |seqidx|
      seq = seqs[seqidx][1]
      seqsize = seq.size
      obs = Array.new(nwords,0)
      options[:wordsize].each{|ws| (0..seqsize-ws).each{|i| wid = wids[seq[i, ws]]; obs[wid] += 1 if not wid.nil?}}
      case options[:scoring_scheme]
      when "bin" then wordobs[seqidx] = obs.map{|x| x > 0 ? 1 : -1}
      when "obs" then wordobs[seqidx] = obs
      end
      
      pbar.inc
    end
    pbar.finish
  end
  
  # sum observed word counts for geneset
  observed = (0..nwords-1).to_a.map{|wid| geneset.inject(0) {|s,sidx| s+wordobs[sidx][wid]}}
  
  begin
    puts "\n>> Estimating expected counts ..."
    pbar = ProgressBar.new("progress",nwords)
    (0..nwords-1).to_a.threach(threads) do |wid|
      # compute expected counts for wid
      exp = (0..bg_samples-1).map{|bgi| seqnbs.inject(0) {|s,nbs| s + wordobs[nbs[bgi]][wid]}}.to_statarray
      emean, esd = exp.mean, exp.stddev
      obs = observed[wid]
      zscores[wid] = [wid,(obs - emean)/esd,obs,"#{emean.to_i}+-#{esd.to_s(1)}"]
      pbar.inc
    end
    pbar.finish
  end
  
else # if bg > 0
  bg = options[:bg]

  puts "\n>> Using sequence background model of order #{bg} ..."

  seqssh = Array.new(nseqs)

  puts "\n>> Shuffling sequences"
  begin
    pbar = ProgressBar.new("progress",nseqs)
    (0..nseqs-1).to_a.threach(threads) do |seqidx| 
      us = UShuffle.new
      us.init_shuffle(seqs[seqidx][1],bg)
      seqssh[seqidx] = (1..bg_samples).to_a.map{|x| seqsh = us.shuffle}
      pbar.inc
    end
    pbar.finish
  end
    
  begin
    puts "\n>> Estimating expected counts ..."
    pbar = ProgressBar.new("progress",nwords)
    wids.to_a.threach(threads) do |word,wid|
      # compute expected counts for wid
      obs = seqs.select{|x| x[1].include?(word)}.size
      exp = bg_samples.times.map{|bgi| nseqs.times.select{|sidx| seqssh[sidx][bgi].include?(word)}.size}.to_statarray
      emean, esd = exp.mean, exp.stddev
      zscores[wid] = [wid,(obs - emean)/esd,obs,"#{emean.to_i}+-#{esd.to_s(1)}"]
      pbar.inc
    end
    pbar.finish
  end
 
end
  
  
###
### FDR
###

# not meaningful here, compute bonferonni corrected p-value from
# normal distribution

###
### Output summarization
###

begin
  spacing = "%-15s"
  wids2 = wids.invert
  zscores = zscores.sort_by{|x| x[1].nan? ? 0.0 : x[1].abs}.reverse
  puts "\nTop #{output_top} words"
  puts ['rank','word','z-score','observed','expected','annotation'].map{|x| sprintf(spacing,x)}.join('')
  zscores[0,output_top].each_with_index do |r,i|
    wd = wids2[r[0]]
    s = [i+1,wd,r[1].to_s(2),r[2],r[3],word_annotation[wd]]
    puts s.map{|x| sprintf(spacing,x)}.join('')
  end

  if options[:dump]
    fname = rankfilename + "." + options[:dump].to_s
    of = File.new(fname,"w")
    of.puts ['rank','word','z-score','observed','expected','annotation'].map{|x| sprintf(spacing,x)}.join('')
    puts "dumping top #{options[:dump]} words in file: #{fname}"
    zscores[0..options[:dump]-1].each_with_index do |r,i|
      wd = wids2[r[0]]
      s = [i+1,wd,r[1].to_s(2),r[2],r[3],word_annotation[wd]]
      of.puts s.map{|x| sprintf(spacing,x)}.join('')
    end
  end
end

