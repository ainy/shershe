#!/usr/bin/pypy
#coding: utf-8
from __future__ import unicode_literals
#import numpy as np

import argparse

parser = argparse.ArgumentParser(description='Align zero.py parser result json files with txt from stdin using Needleman–Wunsch algorithm.')
parser.add_argument('files', metavar='FILE', type=open, nargs='+',
                    help='list of zero.py parser result json files')

parser.add_argument('--debug', dest='debug', action='store_true', default=False,
                    help='debug mode prints unparsed source')

parser.add_argument('--limit', dest='limit', action='store',
                    default=0.25, type=float,
                    help='char error rate limit (default: 0.25)')

args = parser.parse_args()


#Needleman–Wunsch algorithm
insert = 1
delete = 1
match = 0
mismatch = 1

skip_symbol = '#'

def needleman(a,b):
    la = len(a)+1
    lb = len(b)+1
    f = [[0]*lb for i in range(la)]
    for i in range(la):
        f[i][0] = insert*i
    for j in range(lb):
        f[0][j] = delete*j
    for i in range(1,la):
        for j in range(1,lb):
            f[i][j] = min(f[i-1][j-1] + match*(a[i-1]==b[j-1]) + mismatch*(a[i-1]!=b[j-1]), f[i][j-1] + insert, f[i-1][j] + delete)
    
    align_a = ''
    align_b = ''
    i,j=la-1,lb-1
    #print f[i,j]
    while i>0 and j>0:
        if f[i][j] == f[i][j-1] + insert:
            align_a = skip_symbol+align_a
            align_b = b[j-1]+align_b
            j -= 1
        elif f[i][j] == f[i-1][j] + delete:
            align_a = a[i-1]+align_a
            align_b = skip_symbol+align_b
            i -= 1
        else:
            align_a = a[i-1]+align_a
            align_b = b[j-1]+align_b
            i -= 1
            j -= 1
    while i>0:
            align_a = a[i-1]+align_a
            align_b = skip_symbol+align_b
            i -= 1
    while j>0:
            align_a = skip_symbol+align_a
            align_b = b[j-1]+align_b
            j -= 1
    return align_a, align_b, f[la-1][lb-1]/float(len(a))


loct = {1:'первом',2:'втором',3:'третьем',4:'четвёртом',5:'пятом',6:'шестом',7:'седьмом',8:'восьмом',9:'девятом',
        11:'одиннадцатом',12:'двенадцатом',13:'тринадцатом',14:'четырнадцатом',15:'пятнадцатом',16:'шестнадцатом',17:'семнадцатом',18:'восемнадцатом',19:'девятнадцатом',
        10:'десятом',20:'двадцатом',30:'тридцатом',40:'сороковом',50:'пятидесятом',60:'шестидесятом',
        70:'семидесятом',80:'восьмидесятом',90:'девяностом'}
plur = {1:'первых',2:'вторых',3:'третьих',4:'четвёртых',5:'пятых',6:'шестых',7:'седьмых',8:'восьмых',9:'девятых',
        11:'одиннадцатых',12:'двенадцатых',13:'тринадцатых',14:'четырнадцатых',15:'пятнадцатых',16:'шестнадцатых',17:'семнадцатых',18:'восемнадцатых',19:'девятнадцатых',
        10:'десятых',20:'двадцатых',30:'тридцатых',40:'сороковых',50:'пятидесятых',60:'шестидесятых',
        70:'семидесятых',80:'восьмидесятых',90:'девяностых'}
nomn = {1:'первый',2:'второй',3:'третий',4:'четвёртый',5:'пятый',6:'шестой',7:'седьмой',8:'восьмой',9:'девятый',
        11:'одиннадцатый',12:'двенадцатый',13:'тринадцатый',14:'четырнадцатый',15:'пятнадцатый',16:'шестнадцатый',17:'семнадцатый',18:'восемнадцатый',19:'девятнадцатый',
        10:'десятый',20:'двадцатый',30:'тридцатый',40:'сороковой',50:'пятидесятый',60:'шестидесятый',
        70:'семидесятый',80:'восьмидесятый',90:'девяностый'}
accs = {1:'первого',2:'второго',3:'третего',4:'четвёртого',5:'пятого',6:'шестого',7:'седьмого',8:'восьмого',9:'девятого',
        11:'одиннадцатого',12:'двенадцатого',13:'тринадцатого',14:'четырнадцатого',15:'пятнадцатого',16:'шестнадцатого',17:'семнадцатого',18:'восемнадцатого',19:'девятнадцатого',
        10:'десятого',20:'двадцатого',30:'тридцатого',40:'сорокового',50:'пятидесятыго',60:'шестидесятыго',
        70:'семидесятого',80:'восьмидесятого',90:'девяностого'}
quan = {1000:'тысяча',2000:'две тысячи',3000:'три тысячи',4000:'четыре тысячи',5000:'пять тысяч',6000:'шесть тысяч',7000:'семь тысяч',8000:'восемь тысяч',9000:'девять тысяч',
        100:'сто',200:'двести',300:'триста',400:'четыреста',500:'пятьсот',
        600:'шестьсот',700:'семьсот',800:'восемьсот',900:'девятьсот',
        0:'ноль',1:'один',2:'два',3:'три',4:'четыре',5:'пять',6:'шесть',7:'семь',8:'восемь',9:'девять',
        11:'одиннадцать',12:'двенадцать',13:'тринадцать',14:'четырнадцать',15:'пятнадцать',16:'шестнадцать',17:'семнадцать',18:'восемнадцать',19:'девятнадцать',
        10:'десять',20:'двадцать',30:'тридцать',40:'сорок',50:'пятьдесят',60:'шестьдесят',
        70:'семьдесят',80:'восемьдесят',90:'девяносто'}

quan_accs = {1000:'тысячи',2000:'двух тысяч',3000:'трёх тысяч',4000:'четырёх тысячи',5000:'пяти тысяч',6000:'шести тысяч',7000:'семи тысяч',8000:'восми тысяч',9000:'девяти тысяч',
        100:'ста',200:'двухсот',300:'трёхсот',400:'четырёхсот',500:'пятисот',
        600:'шестисот',700:'семисот',800:'восмисот',900:'девятисот',
        0:'ноля',1:'одного',2:'двух',3:'трех',4:'четырех',5:'пяти',6:'шести',7:'семи',8:'восеми',9:'девяти',
        11:'одиннадцати',12:'двенадцати',13:'тринадцати',14:'четырнадцати',15:'пятнадцати',16:'шестнадцати',17:'семнадцати',18:'восемнадцати',19:'девятнадцати',
        10:'десяти',20:'двадцати',30:'тридцати',40:'сорока',50:'пятидесяти',60:'шестидесяти',
        70:'семидесяти',80:'восмидесяти',90:'девяноста'}


orpho = lambda x:[x]

def say_num(num, last=nomn, default=quan):
    say = []
    try:
      if len(num)>=4: say += orpho(default[int(num[-4])*1000])
      if len(num)>=3 and num[-3] != '0': say += orpho(default[int(num[-3])*100])
          
      last0 = last if num[-1]=='0' else default
      if len(num)>=2 and num[-2] == '1':
        say += orpho(last0[int(num[-1])+10])
      else:
        if len(num)>=2 and num[-2] != '0': say += orpho(last0[int(num[-2])*10])
        if len(num)>=1 and num[-1] != '0': say += orpho(last[int(num[-1])])
    except:
      print 'Failed to parse numerical:', num
    return say

import re
def prepare_targets(strkv):
  preps = ['без','из','в','вокруг','из-под','к','меж','над','об','от','перед','под','с','через']
  prepend = ''
  strkv = ' '.join([ st.strip().lower() for st in strkv ]).replace('«','').replace('»','')
  strkv = strkv.replace(' - ',', ').replace(u'\x01',' ').replace(' ',' ').replace('— ',', ').replace('--','')
  strkv = strkv.replace('…','.').replace('?','.').replace(':',',').replace('!','.').replace('...','.').replace(';','.')
  strkv= re.sub(r'\$([0-9,]+) млрд.',r'\1 миллиардов долларов',strkv)
  strkv= re.sub(r'\$([0-9,]+) тысяч',r'\1 тысяч долларов',strkv)
  strkv= re.sub(r'\$([0-9,]+) млн.',r'\1 миллионов долларов',strkv)
  strkv= re.sub(r'\$([0-9,]+)',r'\1 долларов',strkv)
  strkv = strkv.replace(u'трлн.',u'триллионов').replace(u'млрд.',u'миллиардов').replace('млн.','миллионов').replace('тыс.','тысяч').replace('№','номер')
          
  all_targ = []
  ww = strkv.replace(',','').split('.')
  for k in range(len(ww)):
      st = ww[k]
      st= st.strip().replace('*','').replace('"','').replace('%-',' процент').replace('%',' процент').replace('(','').replace(')','')
      st= st.replace('-процентн',' процентн')
      targ = []
      for s in st.split():
          if all([l in 'VIXvix' for l in s]):
              roman = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
              string = s.upper()
              total = 0
              while string:
                  if len(string) == 1 or roman[string[0]] >= roman[string[1]]:
                      total += roman[string[0]]
                      string = string[1:]
                  else:
                      total += roman[string[1]] - roman[string[0]]
                      string = string[2:]
              if prepend: 
                  targ += orpho(prepend)
                  prepend = ''
              targ += say_num(str(total))
              continue
          elif any([l in '0123456789' for l in s]):
              if s[0] == '[': continue #footnote
              if prepend: 
                  targ += orpho(prepend)
                  prepend = ''
              if s.endswith('-м'):
                  targ += say_num(s[:-2], loct)
              elif s.endswith('-и'):
                  targ += say_num(s[:-2], quan_accs, quan_accs)
              elif s.endswith('-е') or s.endswith('-х'):
                  targ += say_num(s[:-2], plur)
              elif s.endswith('-го'):
                  targ += say_num(s[:-3], accs)
              elif '/' in s:
                      a,b = s.split('/')
                      targ += say_num(a, quan)
                      targ += say_num(a, plur)
              elif '-' in s:
                  if s.startswith('-'):
                      targ += orpho('минус')
                      targ += say_num(x)
                  else:
                      for x in s.split('-'): 
                          if any([l not in '0123456789' for l in x]):
                              targ += orpho(x)
                          else:
                              targ += say_num(x)
              else:
                  case = quan
                  if len(s)==4:
                    case = nomn
                  if k+1 < len(ww) and ww[k+1]=='года':
                    case = accs
                  if k+1 < len(ww) and ww[k+1]=='годах':
                    case = plur
                  if k+1 < len(ww) and ww[k+1]=='году':
                    case = loct
                  
                  
                  targ += say_num(s, case)
                  
              continue
          elif all([l in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz' for l in s]):
              #targ += s.lower().rstrip('e')
              continue
          elif not all([l in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz' for l in s]) and \
              any([l in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz' for l in s]):
                  print 'Unusual word must be checked:', s
                  continue
          
          if s:
              if '[' in s: print s
              targ += orpho(s)
              prepend = ''
      if targ: all_targ.append(targ)
  return all_targ

import fcntl, termios, struct, math
rows, columns = struct.unpack('hh', fcntl.ioctl(1, termios.TIOCGWINSZ, '1234'))
def print_needleman(a,b, time):
  cols = int(columns)-1
  for i in range(len(a)/cols):
    er = len([1 for aa,bb in zip(a,b)[i*cols:(i+1)*cols] if aa!=bb])
    print i, '#'*int(math.log(er)), str(time/100/60)+':'+str(time/100%60)
    print a[i*cols:(i+1)*cols]
    print b[i*cols:(i+1)*cols]
    print

#import uniout
import codecs
import json
import sys
#lines  = codecs.open('orig.txt','r','utf-8').readlines()
lines = codecs.getreader('utf-8')(sys.stdin).readlines()
targets = prepare_targets(lines)
targets = '|'.join(map(' '.join, targets))

for f in args.files:
  all_set = []
  data = json.load(f)
  data_r, data_s, data_e = zip(*data)
  data_r = '|'.join(data_r)+'|'
  data_t = []
  print  'processing ', f.name, 'total ',len(data_r),'chars'
  while data_r:
    a,b,ler = needleman(data_r[:5000], targets[:len(data_r[:5000])])
    print '%.5d chars left, сer = %.3f' %(len(data_r), ler)
    b = b.replace('|', ' ')
    if ler > args.limit:
      print 'Error rate limit exceeded, last chunk dumped(json, txt):'
      print_needleman(a,b, time=data_s[len(data_t)])
      
      exit(0)
    aa=a.split('|')
    if len(aa)>2: aa = aa[:-len(aa)/2]
    for a in aa:
      t = b[:len(a)].replace('#','')
      lr = len(a.replace('#',''))+1
      lt = len(t)+1
      
      #print targets[:lt-1]
      #assert t == targets[:lt-1]
      data_r = data_r[lr:]
      targets = targets[lt:]
      b = b[len(a)+1:]
      if lr/(lt or 1) > 2: t=''
      data_t.append(t)
  
  ff = f.name+str('_aligned.json')
  result = filter(lambda x:x[2], zip(data_s[:1]+data_e[:-1], data_e, data_t))
  json.dump(result, codecs.open(ff,'w','utf-8'), ensure_ascii=False)
  if args.debug: print targets[:100]
  """
  while di < len(data):
    t,s,e = data[di]
    di += 1
    targ = ''
    appends = 0
    while len(targ) / float(len(t)) < 0.8:
        if targ: targ += ' '
        targ += ' '.join(targets[ti])
        ti+=1
        appends += 1
        while len(targ) / float(len(t)) > 1.2:
          tt,_,e = data[di]
          di += 1
          t += ' '+tt
          print di,'insert', len(targ) / float(len(t))
          
    a,b,ler = needleman(targ, t)
    
    if len(targ) / float(len(t)) < 1 and ti < len(targets):
      while (len(targ) + len(targets[ti])) / float(len(t)) < 1.2:
        targ2 = targ + ' '.join(targets[ti])
        aa,bb,ler2 = needleman(targ2, t)
        if ler2 < ler:
          print aa
          print bb
          ti += 1
          targ = targ2
          appends += 1
          ler = ler2
          a=aa
          b=bb
          print di,'re-append', ler2, ler
        else: break
    
    if len(targ) / float(len(t)) > 1 and di < len(data):
      while (len(targ)) / float(len(t)+len(data[di][0])) > 0.8:
        t2 = t + data[di][0]
        aa,bb,ler2 = needleman(targ, t2)
        if ler2 < ler:
          e = data[di][2]
          di += 1
          t = t2
          appends += 1
          ler = ler2
          a=aa
          b=bb
          print di,'re-insert', ler2, ler
        else: break
    
    print '%.2d-%.2d %.2f %.2d %.2f'%(di, ti, len(targ) / float(len(t)), appends, ler)
    if (ler > 0.25 or len(targ) / float(len(t))> 1.13) and len(targ) > 10 or ler > 0.5:
      print a
      print b
      raw_input()
    all_set.append((s,e,targ))
    
  json.dump(all_set, open('%.2d.json'%i,'w'))
    
#"""    
