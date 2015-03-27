
major = 0
minor = 1
micro = 0
tag = 'dev'

short = '%d.%d' % (major,minor)
full = '%s%s%s' % (short, '' if tag == '' else '-', tag)
