#!/usr/bin/python

from __future__ import print_function
import sys
from pocketsphinx.pocketsphinx import *
from sphinxbase.sphinxbase import *

from scipy.io import wavfile
from scipy.fftpack import fft

config = Decoder.default_config()
config.set_string('-hmm', '~/cmusphinx-ru-5.2/zero_ru_cont_8k_v3/zero_ru.cd_cont_4000')
config.set_string('-lm', '~/cmusphinx-ru-5.2/zero_ru_cont_8k_v3/ru.lm')
config.set_string('-dict', '~/cmusphinx-ru-5.2/zero_ru_cont_8k_v3/ru.dic')
config.set_string('-logfn', 'sphinx.log')
#config.set_boolean('-mmap', False)
config.set_float('-samprate', 8000.)

config.set_boolean('-remove_noise', False)

print ('init...')
decoder = Decoder(config)
print ('decoding...')
import json
frames = 0
for f in sys.argv[1:]:
  out = []
  stream = open(f, 'rb')
  in_speech_bf = False
  decoder.reinit(config)
  decoder.start_utt()
  while True:
      buf = stream.read(1024)
      if buf:
          decoder.process_raw(buf, False, False)
          if decoder.get_in_speech() != in_speech_bf:
              in_speech_bf = decoder.get_in_speech()
              if not in_speech_bf:
                  decoder.end_utt()
                  out.append ((decoder.hyp().hypstr, list(decoder.seg())[0].start_frame, list(decoder.seg())[-1].end_frame))
                  print (*out[-1])
                  decoder.start_utt()
      else:
          break
  decoder.end_utt()
  
  out.append ((decoder.hyp().hypstr, list(decoder.seg())[0].start_frame, list(decoder.seg())[-1].end_frame))
  print (*out[-1])
                  
  print ('======================',f)
  json.dump(out, open(f+'.json','w'))
