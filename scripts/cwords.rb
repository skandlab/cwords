#!/usr/bin/env jruby

srcdir = File.dirname(__FILE__)
basedir = srcdir + "/../"
libdir = basedir + '/lib/'
$LOAD_PATH << libdir

require 'wordRS-lib.rb'
require 'rubygems'
require 'progressbar'
require 'optparse'
require 'peach'
require 'java'
require libdir + 'ushuffle.jar'
java_import 'UShuffle'

#default options
options = Hash.new
options[:wordsize] = [7]
options[:split_words]=nil
options[:dbdir] = basedir + "db/"
options[:scoring_scheme] = 'pval'
options[:permutations]=50
options[:seqshuffles]=100
options[:rankfile]=nil
options[:seqfile]=nil
options[:report_words]=nil
options[:plot_words]=nil
options[:onlyanno]=nil
options[:dump]=nil
options[:testing]=nil
options[:rank_all]=nil
options[:rank_inverse]=nil
options[:rank_split_median]=nil
options[:rank_abs]=nil
options[:bg]=1 #mononucleotide shuffling
options[:threads]=1

$coptions = OptionParser.new do |opts|
  opts.banner = "Usage: cwords [options]"
  
  # analysis settings
  opts.on("-c", "--scoring_scheme ARG", "scoring scheme") {|o| options[:scoring_scheme] = o} 
  opts.on("-p", "--permutations ARG", "number of list permutations") {|o| options[:permutations] = o.to_i}
  opts.on("-q", "--shuffles ARG", "number of sequence shuffles for sequence bias correction") {|o| options[:seqshuffles] = o.to_i}
  opts.on("-w", "--wordsize ARG", "wordsize") { |o| options[:wordsize] = o.split(",").map{|x| x.to_i}}
  opts.on("-b", "--bg ARG", "background nucleotide model") {|o| options[:bg] = o.to_i}
  opts.on("-t", "--threads ARG", "use multiple threads to parallelize computations") {|o| options[:threads] = o.to_i}
  opts.on(      "--split_words WORDS", "split sequence set based on occurrences of WORDS") {|o| options[:split_words] = o.split(",")}
  opts.on(      "--onlyanno", "only process annotated (i.e. mirbase) words") {|o| options[:onlyanno] = true}
  
  # rank control
  opts.on("-x", "--rank_all", "do not split positive and neg. values") {|o| options[:rank_all] = true}
  opts.on("-m", "--rank_split_median", "split ranked list at median") {|o| options[:rank_split_median] = true}
  opts.on("-i", "--rank_inverse", "inverse all ranked lists") {|o| options[:rank_inverse] = true}
  opts.on("-a", "--rank_abs", "rank by absolute value") {|o| options[:rank_abs] = true}
  
  # files and directories
  opts.on("-r", "--rankfile ARG", "rank file") {|o| options[:rankfile] = o}
  opts.on("-s", "--seqfile ARG", "sequence file") {|o| options[:seqfile] = o}
  opts.on("-d", "--db ARG", "word database") { |o| options[:db] = o}
  
  # output control
  opts.on("-u", "--dump ARG", "dump top words") { |o| options[:dump] = o.to_i}
  opts.on(      "--report_words ARG", "report on words (comma separated)") {|o| options[:report_words] = o.split(',')}
  opts.on(      "--plot_words ARG", "only make plot files for words (comma separated)") {|o| options[:plot_words] = o.split(',')}
  opts.on(      "--testing", "testing mode") {|o| options[:testing] = true}
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
[:rankfile].each{|p| show_help("option '#{p}' mandatory") if options[p].nil?}
show_help("db or seqfile required") if !(options[:db] or options[:seqfile])
show_help("scoring scheme must be one of: obs,bin,pval") if !(['obs','bin','pval'].include?(options[:scoring_scheme]))

testing = options[:testing]

# get filename without directory
rankfilename = File.basename(options[:rankfile])

# hard-coded
output_top = 10

prankdir = basedir + "/db/" + options[:db] + "/" if options[:db]
annofile = basedir + "/resources/" + "word_annotation.tsv" #annotation
tidfile = basedir + "/resources/" + "genemap.tsv"
seqshuffles = 5000 # currently hardcoded for database
sequences = nil
nwords = options[:wordsize].map{|x| 4**x}.to_statarray.sum
bg=options[:bg] # TODO, make option
threads=options[:threads]

###
### Main program
###

puts ">> Parameters"
options.each{|k,v| puts sprintf("%-20s:  %s",k,v) if !v.nil?}

# read in mirbase seed family
word_annotation = Hash.new("") # seq => family
IO.readlines(annofile).each{|l| word_annotation[l.split("\t")[0]] = l.split("\t")[1]}

# read optional sequences
if options[:seqfile]
  puts "\n>> reading sequences ..."
  sequences = Hash.new
  IO.readlines(options[:seqfile],">")[1..-1].each do |entry|
    ls = entry.split("\n").map{|x| x.chomp}
    # hash ensures sequence ids unique
    sequences[ls[0]] = ls[1..-2].join('').downcase.gsub('u','t') # last field is ">"
  end
  seqshuffles = options[:seqshuffles]
end

# initialize word id hash, word sequence => word id (0..nwords-1)
wids = Hash.new
i = 0
options[:wordsize].each{|ws| ['a','g','c','t'].rep_perm(ws) {|seqa| wids[seqa.join('')]=i ; i+=1 }}

###
### ID mapping
###

# pre-computed word database:
#  map ids given in rankfile to internal ids
#  remove rankfile entries with no match to internal id
# sequence file:
#  take intersection of rank and sequence IDs

puts "\n>> Mapping and filtering IDs ..."

all = []
begin
  idmap = Hash.new
  internal_ids = nil

  if sequences
    internal_ids = sequences
  else
    IO.readlines(tidfile).each do |l|
      tid = l.split(" ")[0]
      l.split(" ")[1].split(",").each{|extid| idmap[extid] = tid}
    end
    internal_ids = idmap.invert # allowed internal ids
  end

  allh = Hash.new {|h,k| h[k] = []}
  filtered = 0
  
  IO.readlines(options[:rankfile]).each do |l|
    l = l.split("\t")
    #test if internal id or mapable external id
    tid = (internal_ids.key?(l[0]) ? l[0] :  idmap[l[0]])   
    tid.nil? ? filtered += 1 : allh[tid] << l[1].to_f
  end

  # filter unknown sequences
  sequences.keys.each{|id| sequences.delete(id) if !allh.key?(id)} if sequences
  
  # we currently mean-collapse ids, we could allow mean/min/max collapsing ...
  all = allh.to_a.map{|tid,values| [tid,values.to_statarray.mean]}

  puts "removed #{filtered} invalid transcript ids" if filtered > 0
end

allorder = Hash.new # tid => index in all
all.each_with_index{|x,i| allorder[x[0]] = i}

###
### Word enumeration (optional)
###

wordscores = []
if sequences
  puts "\n>> Enumerating words in sequences"
  wordscores = Array.new(all.size) {Array.new(wids.size,0)} # {Java::short[wids.size].new}
  pbar = ProgressBar.new("progress",sequences.size)
  all.peach(threads) do |seqid,val|
    us = UShuffle.new
    seq=sequences[seqid]
    seqidx=allorder[seqid]
    pbar.inc
    seqsize = seq.size
    observed = Array.new(wids.size,0)
    options[:wordsize].each{|ws| (0..seqsize-ws).each{|i| wid = wids[seq[i, ws]]; observed[wid] += 1 if not wid.nil?}}

    case options[:scoring_scheme]
    when "bin" then wordscores[seqidx] = observed.map{|x| x > 0 ? 1 : -1}
    when "obs" then wordscores[seqidx] = observed
    else
      # pval, compute distribution of expected word occurrences
      us.init_shuffle(seq,bg)
      seqshuffles.times do |si|
        seqsh = us.shuffle
        expected = Array.new(wids.size,0)
        options[:wordsize].each{|ws| (0..seqsize-ws).each{|i| wid = wids[seqsh[i, ws]]; expected[wid] += 1 if !wid.nil?}}
        observed.each_with_index{|x,widx| wordscores[seqidx][widx] =+ 1 if expected[widx]>=x}
      end
    end
  end
  pbar.finish
end

###
### Generate list ranking
###

analyze = []
if options[:rank_split_median]
  # we should perhaps use an :inverse option,
  # reversing the two pos and neg lists
  med = all.map{|x| x[1]}.to_statarray.median
  pos_set = all.select{|x| x[1] > med}.sort{|a,b| b[1] <=> a[1]}
  neg_set = all.select{|x| x[1] <= med}.sort{|a,b| a[1] <=> b[1]}
  analyze = [[pos_set,'med_positive'],[neg_set,'med_negative']]
elsif options[:rank_all] # do not split positive and negative range
  pos_set = all.sort{|a,b| b[1] <=> a[1]}
  neg_set = all.sort{|a,b| a[1] <=> b[1]}
  analyze = [[pos_set,'all_positive'],[neg_set,'all_negative']]
elsif options[:rank_abs] # rank by absolute values
  pos_set = all.map{|x| [x[0],x[1].abs]}.sort{|a,b| b[1] <=> a[1]}
  neg_set = pos_set.reverse 
  analyze = [[pos_set,'abs_positive'],[neg_set,'abs_negative']]
else
  pos_set = all.select{|x| x[1] > 0}.sort{|a,b| b[1] <=> a[1]}
  neg_set = all.select{|x| x[1] < 0}.sort{|a,b| a[1] <=> b[1]}
  analyze = [[pos_set,'positive'],[neg_set,'negative']]
end

# inverse lists
analyze.map!{|set,nm| [set.reverse,nm+".inv"]} if options[:rank_inverse]

# split sequence set when --split option is given
if options[:split_words]
  seqs_with_words = Hash.new

  options[:split_words].each do |split_word|
    begin
      IO.readlines(prankdir + split_word.downcase + ".rnk").each do |x|
        l = x.split("\t")
        seqs_with_words[l[0]] = 1 if l[1].to_i > 0
      end
    rescue
      warn "could not split sequences on word #{split_word}: " + $!
    end
  end
  
  analyze_split = []
  analyze.each do |set,nm|
    analyze_split += set.partition{|x| seqs_with_words.key?(x[0])}.zip([nm+".split+"+options[:split_words].join(","),nm+".split-"+options[:split_words].join(",")])
  end
  analyze = analyze_split
end

###
### Correlation analysis
###

puts "\n>> Analyzing sequence sets: " + analyze.map{|x| x[1]}.join(", ")

analyze.each do |set,nm|
  ngenes = set.size
  puts "\n>> Analyzing #{nm} set ...\nnumber of genes: #{ngenes}"
  next if ngenes == 0
  perms = []
  report = []
  pfdrz = []
  
  franks = Hash.new # tid => index in set
  set.each_with_index{|x,i| franks[x[0]] = i}
  
  puts "permuting #{options[:permutations]} times ...\n"
  options[:permutations].times{|i| perms << (0..set.size-1).to_a.shuffle}
  
  pbar = ProgressBar.new("progress",nwords)
  wids.to_a.sort_by{|x| x[1]}.peach(threads) do |word,wid|
    pbar.inc
    next if options[:onlyanno] and not word_annotation.key?(word) #only process annotated words
    next if options[:plot_words] and !options[:plot_words].include?(word)
    
    plotfile = File.new(rankfilename + ".#{word}.#{nm}.csv","w") if options[:plot_words]
    
    score = Array.new(ngenes) # scores ordered by fold change
    
    if sequences
      score = set.map{|x| wordscores[allorder[x[0]]][wid]}
      score.map!{|x| -Math.log((x+1.0)/(seqshuffles+1))} if options[:scoring_scheme] == 'pval'
    else # use precomputed word database
      wordcounts = IO.readlines(prankdir + word + ".rnk").map{|x| x.split("\t")}.select{|x| franks.key?(x[0])}
      case options[:scoring_scheme]
      when "bin" then wordcounts.each{|id,obs,gte_obs,exp| score[franks[id]] = obs.to_i == 0 ? -1 : 1}
      when "obs" then wordcounts.each{|id,obs,gte_obs,exp| score[franks[id]] = obs.to_f}
      when "pval" then wordcounts.each{|id,obs,gte_obs,exp| score[franks[id]] = -Math.log((gte_obs.to_f+1)/(seqshuffles+1.0))}
      end
    end
    
    smean = score.to_statarray.mean
    maxrs = 0
    leading_edge = 0
    rs = 0 #running sum
    rsa = [0]
    score.each_with_index do |x,i|
      rs += (x-smean)
      rsa << rs
      if rs.abs > maxrs.abs
        maxrs = rs
        leading_edge = i+1
      end
    end
    
    plotfile.puts(([word+".score"] + [0] + score.map{|x| x.to_e(2)}).join(",")) if options[:plot_words]
    plotfile.puts(([word+".rs"] + rsa).join(",")) if options[:plot_words]
    
    # we are only interested in pos. maxrs scores,
    # because we currently analyze up/down regulated seperately
    next if maxrs <= 0
    
    pmaxrs_pos = StatArray.new
    perms.each_with_index do |psa,pidx|
      prs = 0
      prsa = [0]
      pmaxrs = 0
      psa.each do |i|
        prs += score[i]-smean
        prsa << prs
        pmaxrs = prs if prs.abs > pmaxrs.abs
      end
      # the permuted scores are approx. symmetric around 0
      pmaxrs_pos << pmaxrs.abs
      plotfile.puts(([word+".rs."+pidx.to_s] + prsa).join(",")) if options[:plot_words]
    end
    
    pmean = pmaxrs_pos.mean
    pstd = pmaxrs_pos.stddev
    
    #Because the word zscore distr. can be quite different,
    # we compute the deviation from the mean of the absolute dist.
    # The permuted maxRS should be normally distr. (sum of random numbers)
    pfdrz += pmaxrs_pos.map{|x| (x-pmean)/pstd}
    
    #pvalue and fdr statistic for word is also computed based on abs. dist.
    pval = (pmaxrs_pos.select{|x| x>=maxrs}.size+1.0)/(pmaxrs_pos.size+1)
    zsc = (maxrs-pmean)/pstd
    
    plotfile.close if options[:plot_words]
    report << [wid,zsc,pval,nil,leading_edge]
    
  end # wordsize
  pbar.finish
  
  ###
  ### FDR
  ###

  puts "fdr calculation ..."
  fdrrank = pfdrz.map{|x| [x,nil]} # [zscore,word_report_index]
  report.each_with_index{|x,idx| fdrrank << [x[1],idx]}
  fdrrank = fdrrank.sort_by{|x| x[0]}.reverse # sort high zscore to low zscore
  nfp = pfdrz.size.to_f
  ntp = report.size.to_f
  word_fdrrank = Hash.new()
  ifp = 0
  itp = 0
  fdrrank.each do |zsc,idx|
    if idx.nil?
      ifp += 1
    else
      itp += 1
      fpr = ifp/nfp
      tpr = itp/ntp
      report[idx][3] = fpr/tpr
    end
  end
  
  cutoff_fdr = [0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.5]
  puts ""
  puts (["fdr <="] + cutoff_fdr.map{|x| x.to_s(3)} + ["total"]).join("\t")
  puts (["count"] + cutoff_fdr.map{|x| report.select{|y| y[3] <= x}.size} + [report.size]).join("\t")

  ###
  ### Output summarization
  ###
  
  wids2 = wids.invert
  report = report.sort_by{|x| x[1]}.reverse
  puts "\nTop #{output_top} words"
  puts ['rank','word','z-score','p-value','fdr','ledge','annotation'].map{|x| sprintf("%-10s",x)}.join('')
  report[0,output_top].each_with_index do |r,i|
    wd = wids2[r[0]]
    s = [i+1,wd,r[1].to_s(2),r[2].to_e(2),r[3].to_e(2),r[4].to_s,word_annotation[wd]]
    puts s.map{|x| sprintf("%-10s",x)}.join('')
  end

  if options[:report_words]
    puts "......"
    report.each_with_index do |r,i|
      if options[:report_words].include?(r[0]) # and i > output_top
        wd = wids2[r[0]]
        s = [i+1,wd,r[1].to_s(2),r[2].to_e(2),r[3].to_e(2),r[4].to_s,word_annotation[wd]]
        puts s.map{|x| sprintf("%-10s",x)}.join('')
      end
    end
  end

  if options[:dump]
    fname = rankfilename + ".#{nm}." + options[:dump].to_s 
    of = File.new(fname,"w")
    of.puts ['rank','word','z-score','p-value','fdr','ledge','GS size','annotation'].map{|x| sprintf("%-10s",x)}.join('')
    puts "dumping top #{options[:dump]} words in file: #{fname}"
    report[0..options[:dump]-1].each_with_index do |r,i|
      wd = wids2[r[0]]
      s = [i+1,wd,r[1].to_s(2),r[2].to_e(2),r[3].to_e(2),r[4].to_s,word_annotation[wd]]
      of.puts s.map{|x| sprintf("%-10s",x)}.join('')
    end
  end
  
end
