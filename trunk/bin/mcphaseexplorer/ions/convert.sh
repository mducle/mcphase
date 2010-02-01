 #!/bin/bash
 for files in *; do tr \\r \\n < $files > t; mv t $files; done
