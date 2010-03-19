
class Array
  
  # threaded each
  def threach(n = 1, &b)
    return [] if n == 0 or size == 0
    result = Array.new(size)
    self.send(:each,&b) if n == 1
    
    n = [n,size].min
    
    part_size, part_remainder = size/n, size % n
    threads = []

    pstart = 0
    n.times do |pi|
      pend = pstart + part_size - 1
      pend += 1 if pi<part_remainder
      threads << Thread.new(pstart,pend) do |a,b|
        for j in a..b
          yield(slice(j))
        end
      end
      pstart = pend+1
    end
    
    threads.each { |t| t.join }
    self
  end
  
  # unit tests for threach
  #Array.send(:include, Threach)
  #a=(0..4).to_a
  #res = a.map{|x| x*10}
  #a.size.times do |pn|
  #  b=Array.new(a.size)
  #  a.threach(pn+1) {|x| b[x]=x*10}
  #  puts b == beach
  #end

  
  def shuffle()
    arr = self.dup
    arr.size.downto 2 do |j|
      r = Kernel::rand(j)
      arr[j-1], arr[r] = arr[r], arr[j-1]
    end
    arr
  end
  
  #reference: http://blade.nagaokaut.ac.jp/~sinara/ruby/math/combinatorics/array-rep_perm.rb
  def rep_perm(n)
    if n < 0
    elsif n == 0
      yield([])
    else
      rep_perm(n - 1) do |x|
        each do |y|
            yield(x + [y])
          end
      end
    end
  end

  def to_statarray
    StatArray.new(self)
  end

  #return the (sorted) ranks of the elements
  def ranks()
    h = Hash.new {|h,k| h[k]=[]}
    self.sort.each_with_index{|x,idx| h[x] << idx}
    self.map{|x| h[x].first + (h[x].size-1)/2.0}
  end

end

class StatArray < Array

  alias :count size

  def sum
    inject(0) { |sum, x| sum + x }
  end

  def mean
    return 0.0 if self.size == 0
    sum.to_f / self.size
  end
  alias :arithmetic_mean :mean

  def median
    return 0 if self.size == 0
    tmp = sort
    mid = tmp.size / 2
    if (tmp.size % 2) == 0
      (tmp[mid-1] + tmp[mid]).to_f / 2
    else
      tmp[mid]
    end
  end

  # The sum of the squared deviations from the mean.
  def summed_sqdevs
    return 0 if count < 2
    m = mean
    StatArray.new(map { |x| (x - m) ** 2 }).sum
  end

  # Variance of the sample.
  def variance
    # Variance of 0 or 1 elements is 0.0
    return 0.0 if count < 2
    summed_sqdevs / (count - 1)
  end

  # Variance of a population.
  def pvariance
    # Variance of 0 or 1 elements is 0.0
    return 0.0 if count < 2
    summed_sqdevs / count
  end

  # Standard deviation of a sample.
  def stddev
    Math::sqrt(variance)
  end

  # Standard deviation of a population.
  def pstddev
    Math::sqrt(pvariance)
  end

  # Calculates the standard error of this sample.
  def stderr
    return 0.0 if count < 2
    stddev/Math::sqrt(size)
  end

  def gmean
    return 0.0 if self.size == 0
    return nil if self.any?{|x| x < 0 } # not negative
    return 0.0 if self.any?{|x| x == 0 } #includes 0
    return 10**(self.map{|x| log10(x)}.sum/self.size)
  end

  #rank cdf - array is not sorted before cdf calculation
  def rcdf(normalised = 1.0)
    s = sum.to_f
    inject([0.0]) { |c,d| c << c[-1] + normalised.to_f*d.to_f/s }
  end

end


class String
  
  def shuffle
    self.split("").shuffle.join
  end

  def to_f2
    return 1/0.0 if self == "inf"
    return -1/0.0 if self == "-inf"
    return 0/0.0 if self == "NaN"
    self.to_f
  end
  
end


class Float
  
  alias_method :orig_to_s, :to_s

  def to_s(arg = nil)
    if arg.nil?
      #orig_to_s
      sprintf("%f", self)
    else
      sprintf("%.#{arg}f", self)
    end
  end

  def to_e(arg = nil)
    if arg.nil?
      #orig_to_s
      sprintf("%.#{2}e", self)
    else
      sprintf("%.#{arg}e", self)
    end
  end

  def sign
    self < 0 ? -1 : 1
  end
  
end

def log2(number)
  Math.log(number)/Math.log(2)
end

def log10(number)
  Math.log(number)/Math.log(10)
end

