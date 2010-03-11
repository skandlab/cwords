#!/usr/bin/ruby

###
### Given a miRS output word list file and a sequence,
###  the words are aligned to input sequence (RNA complementary alignment, differentiating between mismatches and bulges)
###

srcdir = File.dirname(__FILE__)
$LOAD_PATH << srcdir + '/../lib/'

require 'wordRS-lib.rb'
require 'progressbar'
require 'optparse'
require 'pp'

#default options
options = Hash.new
options[:wordfile]=nil
options[:sep]=" "
options[:overlap]=3
options[:seedoverlap]=nil
options[:testing]=nil
options[:fdr]=nil
options[:seed]=20
options[:top]=nil
options[:keep_lc] = nil # filter low complexity words
# we could estimate significance of cluster size based on shuffles ...
options[:shuffles]=100 

$coptions = OptionParser.new do |opts|
  opts.on("-w", "--wordfile ARG", "word rank file") {|o| options[:wordfile] = o}
  opts.on("-s", "--sep ARG", "separator") {|o| options[:sep] = o}
  opts.on("-k", "--keep_lc", "keep low complexity words") {|o| options[:keep_lc] = o.to_i}
  opts.on("-o", "--overlap ARG", "overlap") {|o| options[:overlap] = o.to_i}
  opts.on("-v", "--seedoverlap ARG", "seed overlap") {|o| options[:seedoverlap] = o.to_i}
  opts.on("-m", "--matchscore ARG", "scoring scheme") {|o| options[:matchscore] = o.to_f}
  opts.on("-t", "--top ARG", "limit words to scan") { |o| options[:top] = o.to_i}
  opts.on("-i", "--seed ARG", "seed number clusters") { |o| options[:seed] = o.to_i}
  opts.on("-f", "--fdr ARG", "limit words to scan") { |o| options[:fdr] = o.to_f}
  opts.on("-u", "--shuffles ARG", "background estimate based on shuffle sequences") { |o| options[:shuffles] = o.to_i}
end

def show_help(msg="", code=0, io=STDOUT)
  io.puts "#{msg}\n#{$coptions}"
  exit(code)
end

$coptions.parse!(ARGV)
#mandatory parameters
[:wordfile].each{ |p| show_help("option '#{p}' mandatory") if options[p].nil?}

options[:seedoverlap] = options[:overlap] if !options[:seedoverlap]
options[:top] = options[:seed] if !options[:top]

# read in mirbase seed family
annofile = srcdir + "/../resources/" + "word_annotation.tsv" #annotation
@word_annotation = Hash.new("") # seq => family
IO.readlines(annofile).each{|l| @word_annotation[l.split("\t")[0]] = l.split("\t")[1]}

###
### Sub
###

class String
  
  def editdist(b)
    a = self
    return 0 if !a || !b || a == b
    return (a.length - b.length).abs if a.length == 0 || b.length == 0
    m = [[0]]
    1.upto(a.length) {  |i| m[i] = [i] }
    1.upto(b.length) {  |j| m[0][j] = j }
    1.upto(a.length) do |i|
      1.upto(b.length) do |j|
        m[i][j] =
          [ m[i-1][j-1] + (a[i-1] == b[j-1] ? 0 : 1),
            m[i-1][j] + 1,
            m[i][j-1] + 1                             ].min
      end
    end
    m[a.length][b.length]
  end

  def local_align_score(b)
    a = self
    return 0 if !a || !b || a == b
    return (a.length - b.length).abs if a.length == 0 || b.length == 0
    match = 2
    miss = -1
    gap = -1
    m = [[0]]
    1.upto(a.length) {  |i| m[i] = [0] }
    1.upto(b.length) {  |j| m[0][j] = 0 }
    1.upto(a.length) do |i|
      1.upto(b.length) do |j|
        m[i][j] =
          [ m[i-1][j-1] + (a[i-1] == b[j-1] ? match : miss),
            m[i-1][j] + gap,
            m[i][j-1] + gap,0                             ].max
      end
    end
    (m.last+m.map{|x| x.last}).max
  end

  # gapped align of complementary nucleotides, with wobbles
  def local_align_gap_score(b)
    a=self
    return 0 if !a || !b || a == b
    return (a.length - b.length).abs if a.length == 0 || b.length == 0
    a = a.split('')
    b = b.split('')

  # score cutoff 8.5
  # at least 5 matches, max 1 go or mismatch
  # wobbles and neutral
  # allowed
  # 6m + 1go/mm : 11.0
  # 5m + 2wo : 10
  # 5m + 1go/mm + 1wo : 9
  # 5m + 1go + 1ge : 8.5
  # not allowed
  # 4m + 3wo : 8 # no
  # 5m + 2mm/go : 8

    match = 2.0
    miss = -1.0
    wobble = 0.0
    gap_open = -1.0
    gap_ext = -0.5

    score = Hash.new(miss)
    ['at','ta','gc','cg'].each {|x| score[x] = match}
    ['gt','tg'].each {|x| score[x] = wobble}

    m = [[0]]
    #g1/2 - gap extension matrix for two seqs
    g1 = [[gap_open]]
    1.upto(a.length) {  |i| m[i] = [0]; g1[i] = [gap_open]} #g2[i] = gap_open}
    1.upto(b.length) {  |j| m[0][j] = 0; g1[0][j] = gap_open}# g2[0][j] = gap_open}
    g2 = g1.clone

    1.upto(a.length) do |i|
    1.upto(b.length) do |j|
      scores = [ m[i-1][j-1] + score[a[i-1] + b[j-1]], #match/mismatch
                 m[i-1][j] + g1[i-1][j],  # gap in seq1
                 m[i][j-1] + g2[i][j-1],  # gap in seq2
                 0] # start new local alignment

        m[i][j] = scores.max
        g1[i][j] = (scores.index(m[i][j]) == 1) ? gap_ext : gap_open # gap in seq1
        g2[i][j] = (scores.index(m[i][j]) == 2) ? gap_ext : gap_open # gap in seq2
      end
    end
    #puts ([""]+b).map{|x| x.center(4)}.join(' ')
    #puts m.map{|x| x.map{|y| y.to_f.to_s(1).center(4)}.join(" ")}
    (m.last+m.map{|x| x.last}).max
  end
  
  # returns [size of overlap,rel. pos of s2 in s1]
  def overlap(s2,require_ovl)
    s1 = self
    s1=s1.split('')
    s2=s2.split('')
    ovl = require_ovl || 5
    return [] if (s1.size < ovl or s2.size < ovl)
    return [0] if s1 == s2
    
    # first, check if one word is contained in the other ...
    # if s1.size <= s2.size
    (s2.size-s1.size).times{|i| return [-i] if s2[i,s1.size] == s1}
    (s1.size-s2.size).times{|i| return [i] if s1[i,s2.size] == s2}
        
    # then check for overlap
    (ovl..[s1.size,s2.size].min).to_a.reverse.each do |len|
      ovls = []
      ovls << -(s2.size-len) if s1[0,len] == s2[s2.size-len..-1]
      ovls << s1.size-len if s2[0,len] == s1[s1.size-len..-1]
      return ovls if !ovls.empty?
    end
    
    return [] # no overlap
  end
  
  #require '~/miwork/words/src/lib/mirs-lib.rb'
  def editaln(b,edt)
    a = self
    #puts a
    #puts b
    return nil if !a || !b
    return nil if a.length == 0 || b.length == 0
    m = [[0]]
    bt = [[0]] # backtrack
    1.upto(a.length) {  |i| m[i] = [i]; bt[i] = [0];  }
    1.upto(b.length) {  |j| m[0][j] = j; bt[0][j] = 0;}
    1.upto(a.length) do |i|
      1.upto(b.length) do |j|
        scores = 
          [ [m[i-1][j-1] + (a[i-1] == b[j-1] ? 0 : 1),0],
            [m[i-1][j] + 1,1],
            [m[i][j-1] + 1,2]      ]
        
        m[i][j] = scores.min.first
        bt[i][j] = scores.min.last
        
      end
    end
    editdist = m[a.length][b.length]
    return [] if editdist > edt
    #find alignment
    if bt[a.length][b.length] == 0
      # last pos aligns
      return [a.length-b.length]
    elsif bt[a.length][b.length] == 1
      # aligns from above, find first non-1
      bt.map{|x| x.last}.reverse.each_with_index do |bts,idx|
        return [a.length-idx-b.length] if bts != 1
      end   
    else
      # aligns from left, find first non-2
      bt[a.length].reverse.each_with_index do |bts,idx|
        return [a.length+idx-b.length] if bts != 2
      end      
    end
        
    puts "error"
    pp m
    pp bt

  end
  
end



def align_to_cluster(cluster,s,ovl)
  align = []
  return [] if cluster.key?(s) # s already exist in cluster
  match_to_cluster = 0
  cluster.to_a.each do |w,wis|
    wol = w.overlap(s,ovl)
    #wol = w.editaln(s,ovl) # test editdist
    if !wol.empty?
      match_to_cluster += 1
      wol.each do |si|
        wis.each do |wi|
          #compute overlap length
          olength = (si < 0) ? [s.size+si,w.size].min : [s.size,w.size-si].min
          # we require the overlap should contain at least two
          # different nucleotides
          #puts w if s == "tgttt" and olength >= 3
          next if (si<0) and w[0,olength].split('').uniq.size == 1
          next if (si>=0) and s[0,olength].split('').uniq.size == 1
          #adjust si relative to first cluster word
          align << [s,wi+si,olength]
        end
      end
    end
  end
  # return [] if overlap is not greater than ovl in at least half of
  # the cluster members
  #return [] if cluster.size > match_to_cluster
  
  # return alignments with greatest overlap
  align = align.compact.uniq
  maxovl = align.map{|x| x.last}.max
  return align.select{|x| x.last == maxovl}
end

# find cluster(s) with highest overlap
def align_to_clusters(clusters,s,ovl)

  alns = []
  clusters.each do |cl|
    aln = align_to_cluster(cl,s,ovl)
    alns << [aln,cl] if !aln.empty?
  end
  #pp alns if s == "tgttt"

  # select alignments with greates overlap
  maxovl = alns.map{|aln,cl| aln.first.last}.max
  #pp maxovl if s == "cccgttt"
  #pp alns.select{|x| x.first.last == maxovl} if s == "cccgttt"
  return alns.select{|aln,cl| aln.first.last == maxovl}
  
end

def add_to_clusters(alignments)
  #add.compact.uniq.each{|x,y| cl[x] << y}
  alignments.each do |alns,cl|
    #pp alns
    alns.each{|aln| cl[aln[0]] << aln[1]}
  end
end

def print_cluster(cluster,ranks)
  cl = []
  cluster.to_a.each{|a,b| b.each{|c| cl << [a,c]}}
  cl = cl.sort_by{|x| x[1]}
  imin = cl.first[1]

  # preprocess clusters for consensus alignment
  cons = Hash.new {|h,k| h[k] = Hash.new(0.0)}
  cl.each do |w,i|
    w.split('').each_with_index do |nt,idx|
      cons[i+idx][nt] +=1
    end
  end

  def cons_word(consh,word,i)
    cw = ""
    word.split('').each_with_index do |nt,idx|
      frac = consh[i+idx][nt]/consh[i+idx].values.to_statarray.sum
      cw += (frac >= 0.5 and consh[i+idx].values.to_statarray.sum > 1) ? nt.upcase : nt
    end
    return cw
  end
  
  cl.each do |w,i|
    puts ranks[w].to_s.rjust(5)+": "+(" "*(i-imin)+cons_word(cons,w,i).tr('tT','uU')).ljust(20)+@word_annotation[w].chomp.ljust(20)
  end
end

###
### Main
###

aw = Hash.new()
if options[:fdr]
  IO.readlines(options[:wordfile]).select{|x| x.split(' ')[5].to_f < options[:fdr]}[0,options[:top]+1][1..-1].each_with_index do |wline,idx|
    aw[idx+1] = wline.split(options[:sep])[1]
  end
else
  IO.readlines(options[:wordfile])[0,options[:top]+1][1..-1].each_with_index do |wline,idx|
    aw[idx+1] = wline.split(options[:sep])[1]
  end
end

awords = aw.to_a.sort.map{|x| x[1]}
awords_seed = awords[0,options[:seed]]
awords_rest = awords[options[:seed]..-1]

if !options[:keep_lc]
  awords_seed = awords_seed.select{|x| x.split('').uniq.size > 1}
  awords_rest = awords_rest.select{|x| x.split('').uniq.size > 1}
  puts "removed #{aw.size-awords_seed.size-awords_rest.size} low complexity words"
end

# seedialize clusters
clusters = []
awords_seed.each{|x| h = Hash.new{|j,k| j[k]=Array.new}; h[x]=[0]; clusters << h}
#pp clusters

puts "step 1, seedializing clusters from #{options[:seed]} words"
pbar = ProgressBar.new("running",options[:seed])
options[:seed].size.times do  
  awords_seed.each do |word|
     add_to_clusters(align_to_clusters(clusters,word,options[:seedoverlap]))
  end
  pbar.inc
end
pbar.finish

# remove duplicate clusters
# {'word' => [pos in alignment,array, usually one element]}

clusters.each{|x| clusters.delete_if{|y| x.object_id != y.object_id and x.key?(y.to_a.first.first)}}

# pp clusters
# merge clusters ...

puts "step 2, processing #{awords.size} words"
#pbar = ProgressBar.new("running",awords.size)
#awords[options[:seed]..-1].size.times do  

# map extra words to cluster(s) with greatest overlap
# dont extend extra words
if options[:top] > options[:seed]
  add = []
  awords_rest.each do |word|
    #    add += align_to_cluster(cl,word,options[:overlap])
    add << align_to_clusters(clusters,word,options[:overlap])
  end
  add.each{|a| add_to_clusters(a)}
end

#require 'pp'
#pp clusters.select{|x| x.size > 2}
#pp clusters
wa = aw.invert
resc = clusters.select{|x| x.size >= 3}
resc.each{|cl| print_cluster(cl,wa);puts "\n"}

puts "Found #{resc.size} word clusters."

