# Ubuntu 
FROM ubuntu:24.04

# Add procps to allow the use of ps command in the container
RUN apt-get update && apt-get install -y perl procps && rm -rf /var/lib/apt/lists/*

# Add the Phobius package to the container. add command will extract the contents of the tarball into /opt/phobius directory.
ADD phobius101_linux.tgz /opt/

# ensure the phobius.pl script and all files in the directory are executable for all users and change to the 64bit version of the decodeanhmm binary in the phobius.pl script.
RUN chown -R root:root /opt/phobius && \
    chmod -R a+rx /opt/phobius && \
    sed -i 's/"$PHOBIUS_DIR\/decodeanhmm"/"$PHOBIUS_DIR\/decodeanhmm.64bit"/' /opt/phobius/phobius.pl

# add the phobius directory to the PATH environment variable so that the phobius.pl script can be run from anywhere in the container.
ENV PATH="/opt/phobius:${PATH}"