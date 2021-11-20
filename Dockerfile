FROM wbarshop/msptodb_base
MAINTAINER Direct All Questions Directly To Devin Schweppe "dkschwep@uw.edu"

#Stage Shiny files, etc.
RUN mkdir /root/Repos/
COPY . /root/Repos/MSPtoDB/

EXPOSE 3838

CMD ["Rscript", "/root/Repos/MSPtoDB/shiny/app.R"]
