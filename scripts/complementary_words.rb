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

#default options
options = Hash.new
options[:wordfile]=nil
options[:seq]=nil
options[:matchscore]=7
options[:testing]=nil
options[:fdr]=nil
options[:top]=250
options[:shuffles]=100

$coptions = OptionParser.new do |opts|
  opts.on("-w", "--wordfile ARG", "word rank file") {|o| options[:wordfile] = o}
  opts.on("-s", "--seq ARG", "sequence") {|o| options[:seq] = o}
  opts.on("-m", "--matchscore ARG", "scoring scheme") {|o| options[:matchscore] = o.to_f}
  opts.on("-t", "--top ARG", "limit words to scan") { |o| options[:top] = o.to_i}
  opts.on("-f", "--fdr ARG", "limit words to scan") { |o| options[:fdr] = o.to_f}
  opts.on("-u", "--shuffles ARG", "background estimate based on shuffle sequences") { |o| options[:shuffles] = o.to_i}
end

def show_help(msg="", code=0, io=STDOUT)
  io.puts "#{msg}\n#{$coptions}"
  exit(code)
end

$coptions.parse!(ARGV)
#mandatory parameters
[:wordfile,:seq].each{ |p| show_help("option '#{p}' mandatory") if options[p].nil?}
seq = options[:seq]

###
### Sub
###

def local_align(s1,s2)
  return 0 if !s1 || !s2 || s1 == s2
  return (s1.length - s2.length).abs if s1.length == 0 || s2.length == 0
  a = s1.split('')
  b = s2.split('')

  # a, word, short seq, target sequence : gaps (bulge in miR) expensive, but also because we have additional match options
  # b, miRNA, longer seq,
  # we require full alignment without penalty for overhangs

  # score cutoff 8.5
  # at least 5 matches, max 1 go or mismatch
  # wobbles and neutral
  # gm : gap mir, gt = gap target
  # allowed
  # 6m + 1gmo/mm : 11.0
  # 5m + 2wo : 10
  # 5m + 1gmo/mm + 1wo : 9
  # 5m + 1gmo + 1ge : 8.5
  # not allowed
  # 4m + 3wo : 8 # no
  # 5m + 2mm/go : 8

  match = 2.0
  miss = -1.0
  wobble = 0.0
  g1_open = -1.0 #
  g2_open = -3.0 # 7m + 1g2 = 14-3 = 11 = 6m + 1g1 = 12 - 1
  gap_ext = -0.5

  score = Hash.new(miss)
  ['at','ta','gc','cg'].each {|x| score[x] = match}
  ['gt','tg'].each {|x| score[x] = wobble}

  m = [[0]] # score matrix
  bm = [["1,1:"]] # backtrack matrix
  g1 = [[g1_open]]
  g2 = [[g2_open]]
  1.upto(a.length) {  |i| m[i] = [0]; g1[i] = [g1_open]; g2[i] = [g2_open]; bm[i]=["#{i+1},1:"]}
  1.upto(b.length) {  |j| m[0][j] = 0; g1[0][j] = g1_open; g2[0][j] = g2_open; bm[0][j]="1,#{j+1}:"}
  #g1/2 - gap extension matrix for two seqs

  1.upto(a.length) do |i|
    1.upto(b.length) do |j|
      scores = [ m[i-1][j-1] + score[a[i-1] + b[j-1]], #match/mismatch
                 m[i-1][j] + g1[i-1][j],  # gap in seq1
                 m[i][j-1] + g2[i][j-1]]  # gap in seq2

      m[i][j] = scores.max
      case scores.index(m[i][j])
      when 0
        g1[i][j] = g1_open
        g2[i][j] = g2_open
        case score[a[i-1] + b[j-1]]
        when match then bm[i][j] = bm[i-1][j-1] + "*"
        when miss then bm[i][j] = bm[i-1][j-1] + "!"
        when wobble then bm[i][j] = bm[i-1][j-1] + "w"
        end
      when 1
        g1[i][j] = gap_ext
        g2[i][j] = g2_open
        bm[i][j] = bm[i-1][j] + "-"
      when 2
        g2[i][j] = gap_ext
        g1[i][j] = g1_open
        bm[i][j] = bm[i][j-1] + "_"
      end

    end
  end

  #  puts ([""]+b).map{|x| x.center(4)}.join(' ')
  #  puts m.map{|x| x.map{|y| y.to_f.to_s(1).center(4)}.join(" ")}
  #(m.last+m.map{|x| x.last}).max
  (m.last+m.map{|x| x.last}).zip(bm.last+bm.map{|x| x.last}).sort_by{|x| x[0]}.last
end

def update_mirpairing(mphash,aln)
  pos,alnstr = aln.split(":")
  start_word,start_mir = pos.split(",")
  alnarr = alnstr.split('')
  mirpos = start_mir.to_i
  while !alnarr.empty?
    alnchar = alnarr.shift
    mphash[alnchar][mirpos-1] += 1
    mirpos += 1 if alnchar != "-" #next mirpos, unless gap in miR
  end
end

###
### Main
###

# this is the string (the miRNA target sequence) we want to match
seqr = seq.downcase.reverse
seqc = seq.downcase.tr("agctu","tcgaa") # control - miR*
shuffled_seqs = (1..options[:shuffles]).to_a.map{|i| seqr.shuffle}
matches_shuffles = Array.new(options[:shuffles],0)
random_seqs = (1..options[:shuffles]).to_a.map{|i| ("agct".split('')*100).shuffle[0,seq.size].join('')}
matches_random = Array.new(options[:shuffles],0)

ofname = "#{options[:wordfile]}.#{seq}"
of = File.new(ofname,"w")
of.puts ">revcompl_sequence\n#{seqc.reverse}"

matches_seq = 0 # matches to the original string

aln_flank = 0 # how should we add flank info ??
mirpairing = Hash.new()
mirpairing['*'] = Array.new(seq.size+aln_flank,0) # match
mirpairing['!'] = Array.new(seq.size+aln_flank,0) # mismatch
mirpairing['w'] = Array.new(seq.size+aln_flank,0) # wobble
mirpairing['-'] = Array.new(seq.size+aln_flank,0) # gap in miR, bulge in target - hard to visualize
mirpairing['_'] = Array.new(seq.size+aln_flank,0) # gap in target, bulge in miR

awords = []
if options[:fdr]
  awords = IO.readlines(options[:wordfile]).select{|x| x.split(' ')[5].to_f < options[:fdr]}[0,options[:top]+1]
else
  awords = IO.readlines(options[:wordfile])[0,options[:top]+1]
end
awords.shift # remove header

puts "Processing #{awords.size} words"
pbar = ProgressBar.new("running",awords.size)
awords.each_with_index do |wline,idx|
  rank = idx+1
  word = wline.split(" ")[1]
  shuffled_seqs.each_with_index{|s,idx| matches_shuffles[idx] += 1 if local_align(word,s).first >= options[:matchscore]}
  random_seqs.each_with_index{|s,idx| matches_random[idx] += 1 if local_align(word,s).first >= options[:matchscore]}
  sc,aln = local_align(word,seqr)
  if sc >= options[:matchscore]
    matches_seq += 1
    update_mirpairing(mirpairing,aln)
    of.puts ">#{rank}[#{sc.to_s(1)}=#{aln}]\n#{word}"
  end
  pbar.inc
end
pbar.finish
of.close

shuffled_mean = matches_shuffles.to_statarray.mean
shuffled_stdd = matches_shuffles.to_statarray.stddev
random_mean = matches_random.to_statarray.mean
random_stdd = matches_random.to_statarray.stddev

puts "words reverse-complementary to sequence: #{matches_seq}"
puts "words reverse-complementary to shuffled : #{shuffled_mean.to_s(1)}+-#{shuffled_stdd.to_s(1)}"
puts "words reverse-complementary to random   : #{random_mean.to_s(1)}+-#{random_stdd.to_s(1)}"
puts "similar words in file: #{ofname}"

puts (["3'"] + seqc.reverse.split('')).map{|x| x.center(4)}.join('')
mirpairing.each do |alnchar,mirpos|
  alnstr = alnchar.center(4)
  alnstr += mirpos.map{|x| x.to_s.center(4)}.join('')
  puts alnstr
end
