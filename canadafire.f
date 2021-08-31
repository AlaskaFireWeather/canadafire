C NCLFORTSTART
        subroutine FIRECODES (tin,hin,win,rin,buiout,fmout,
     &                       isiout,fwiout,dmcout,dcout)
C        integer tlen
        real tin(183),hin(183)
        real win(183),rin(183)
        real fmout(183),isiout(183)
        real fwiout(183),dmcout(183)
        real buiout(183),dcout(183)
C        integer indate,yr,mon,day,hr
C        character*10 timez
C NCLEND
C        write(timez,'(i10)') indate
C        read(timez,'(i4,3i2)') yr,mon,day,hr

C        call main1io(t2,precp,q2,press,yr,mon,day,hr,tlen,bui)
        call canadafire(tin,hin,win,rin,buiout,fmout,isiout,
     &           fwiout,dmcout,dcout)
        
        return 
        end
c***************************************************************************************
        subroutine canadafire(tin,hin,win,rin,buiout,fmout,isiout,
     &        fwiout,dmcout,dcout)
c input an array of tin, hin, win, and rin from one grid point, from April 1 to Sept 30, julian day 91 (april 1) to julian day 273
c tin is in C, hin is %, win km/hr, and rin in mm
        implicit none
c
c Standard 1984 version of Canadian Forest Fire Weather Index System 
c as of 30 August 1984
c recoded by U Bhatt 2/24/2016 updated FORTRAN
c
c Written in FORTRAN 77 for DEC-PDP-11/44 at PNFI
c
c Reads weather data in either english or metric units and prints out in metric only
c 
c codes and indices output to one decimal place.
c daily severity rating output to two decimal places.
c
c LMON = Length of Months
c EL   = DMC Day Length Factors
c FL   = DC Day Length Factors

c declarations
c array of data from the main program
        real tin(91:273), hin(91:273),win(91:273),rin(91:273)
        real fmout(91:273), isiout(91:273),fwiout(91:273),dmcout(91:273)
        real buiout(91:273),dcout(91:273)
c 
        real le(12),lf(12)
        real undef
        integer lmon(12)
        real t, w, r, to,ho
        real h
        integer im, iday,iswitch,jday,nday
c fine fuel moisture code variables
         real mo, rf, mr, ed,ew, ko,kd,kl,kw,m,ffm
c build up index 
        real ffmc0,dmc0,dc0
        real b,bui,u,v,val1,d,dc,dmc
        real dr,dsr,fd,ff,fo,fw,fwi,isi,k
        real p,pd,po,pold,pr,qo,qr,rd,re,s
        character*50 junkie


c data statement
       data lmon/ 31,28,21,30,31,30,31,31,30,31,30,31/
cc ## this is twice the value in instructions
       data le/6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0/
       data lf/-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5.0,2.4,0.4,-1.6,-1.6/

c&& enter in the inital values here
c 3 fuel moisture codes, initialize
       ffmc0 =85.5
       dmc0=6.0
       dc0=15.0
       
       

c keep track of julian day
            jday=90
c start on april 1 and go to september 30 
       do 290 im=4,9  
         do 290 iday=1,lmon(im)
            jday = jday +1
c move the data from the array to the variable for the calculation on a given day
            t = tin(jday)
            h = hin(jday)
            w = win(jday)
            r = rin(jday)

C           write (6,*) im, lmon(im),iday,jday
cc write in a days of data t, h,w,r
C         write(6,*) jday, nday, t, h, w, r
   
 
cc calculate output file that contains
c data t h w r  and ffmc dmc dc isi bui fwi dsr

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c fine fuel moisture code 
cc EQ 1
         mo = (147.2*(101.-ffmc0))/(59.5+ffmc0)   
c calculate rf
cc EQ 2
         rf=r
         if(rf.gt.0.5) then
            rf=r-0.5    

c If r<0.5 then you can not use the following equations 3a and 3b
c then the rainfall routine must be omitted
cc Condition 1 
c calculate wmo (fine fuel moisture code from previous day, mo) 
            if(mo.le.150.)then
cc EQ 3a
                mr= mo + 42.5*rf*exp(-100./(251.-mo))*
     *           (1.-exp(-6.93/rf))
            else
cc EQ 3b
                mr=mo+42.5*rf*exp(-100./(251.-mo))*(1.-exp(-6.93/rf))+
     *        (0.0015*(mo-150.)**2)*rf**0.5
            endif
cc Condition 3
            if(mr.gt.250.) mr=250.
              mo=mr
         endif 
c calculate ed
cc EQ 4 
         ed=0.942*(h**0.679)+11.*exp((h-100.)/10.)+0.18*(21.1-t)
     *     *(1.-1./exp(0.115*h))

c       iswitch=0

       if(mo.gt.ed)then
cc EQ 6a
       ko = 0.424*(1.-(h/100)**1.7)+(0.0694*(w**.5))*
     *       (1.-(h/100.)**8)
cc EQ 6b
       kd = ko* 0.581*exp(0.0365*t)
cc EQ 8 
       m = ed+(mo-ed)*10**(-1.*kd)
c       iswitch=1
       else

cc EQ 5 
        ew=0.618*(h**0.753)+10.*exp((h-100.)/10.)+0.18*(21.1-t)
     *       *(1.-exp(-0.115*h))
         if(mo.lt.ew)then
cc EQ 7a
           kl=0.424*(1.-((100.-h)/100)**1.7)+(0.0694*(w**.5))*
     *   (1.-((100.-h)/100.)**8)
cc EQ 7b
           kw =kl*(0.581*exp(0.0365*t))
cc EQ 9 
           m=ew-(ew-mo)*10.**(-1.*kw)
         else
           m=mo
         endif
c      iswitch=1
       endif

cc direction 8
c        if(iswitch.eq.0) m=mo
cc Condition 2 
c        if(m.gt.250.) m=250.

c calculate fine fuel moisture code
cc EQ 10
       ffm=(59.5*(250.-m))/(147.2+m)
cc restrictions on ffm
       if(ffm.gt.101.)ffm=101.
       if(ffm.le.0.)ffm=0.  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c
c Duff moisture code
c
c yesterdays code becomes dmc0
c direction 1
           po=dmc0

c condition 1 r must be above zero to calculate these quantities. 
         if(r.le.1.5)then
             pr =po
         else
cc EQ 11
           re=0.92*r-1.27
cc EQ 12
           mo= 20.0 + 280.0/ exp(0.023*po)
c----------------------------------------------------------
c calculate b depending on po value
           if(po.le.33.)then
cc EQ 13a
             b=100./(0.5+0.3*po)
           else
cc EQ 13c
             b=6.2*log(po)-17.2 
             if(po-65.0.le.0) b=14.-1.3*log(po) 
           endif
cc EQ 13b
c----------------------------------------------------------
cc EQ 14
          mr=mo+(1000.*re)/(48.77+b*re)
cc EQ 15
          pr=43.43*(5.6348-log(mr-20.))

         endif
         if(t.ge.-1.1)then
          k=1.894*(t+1.1)*(100.-h)*(le(im)*0.0001)
         else
          k=0.0
         endif
c condition 2 that pr can not be less than zero
         if(pr.lt.0.)pr=0.
           dmc = pr+k
c condition 3 t must be greater than -1.1       
         if(dmc.le.0)dmc=0.

 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c drought code
c
cc condition 1
             if(r.gt.2.8) then     
cc EQ 18  
                rd=0.83*r-1.27
cc EQ 19  
                qo=800.*exp(-dc0/400.)
cc EQ 20  
                qr=qo+3.937*rd             
cc EQ 21  
                dr = 400.*log(800./qr)
                 if(dr.gt.0.0)then
                    dc0=dr
                 else
                    dc0=0.0
                 endif
             endif

cc condition 3
             if(t.lt.-2.8)then
              v=lf(im)
             else
cc EQ 22
              v = 0.36*(t+2.8) + lf(im)
cc condition 4
             endif
             if(v.le.0)v=0.
cc EQ 23
              d = dc0 + 0.5*v        

cc set drought code
             dc=d


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c
c initial spread index, buildup index, fire weather index
c initial spread index
cc EQ 24
              
              fw = exp(0.05039*w) 
cc EQ 24
              ff = 91.9*exp(-0.1386*m)*(1.+(m**5.31)/4.93e07)
cc EQ 26
cc initial spread index
              isi = 0.208*fw*ff

c buildup index
cc EQ 27a
             if(dmc.le.0.4*dc)then
                 u = 0.8*dmc*dc/(dmc+0.4*dc) 
             else
cc EQ 27b
                u = dmc - (1. -0.8*dc/(dmc+0.4*dc))*
     *         (0.92 + (0.0114*dmc)**1.7)
             endif
c buildup index
              if (u.lt.0)u=0.0

c fire weather index
cc EQ 28a
             if(u.le.80.) then
              fd=0.626*u**0.809 + 2.
             else
cc EQ 28b
              fd=1000./(25. + 108.64*exp(-0.023*u))
             endif
cc EQ 29 
             b= 0.1*isi*fd
cc EQ 30a
             if(b.gt.1.)then
               s=(exp(2.72*(0.434*log(b))**0.647))
             else
               s=b
             endif
cc EQ 30b
c todays fire weather index
             fwi=s
           
c daily severity rating
             dsr=0.0272*fwi**1.77

c housekeeping items
c set todays values to the yesterdays values before going on
                 ffmc0=ffm
                 dmc0=dmc
                 dc0=dc

c move the bui array item into the array for the day
                 buiout(jday)=u
                 fmout(jday)=ffm
                 isiout(jday)=isi
                 fwiout(jday)=s
                 dmcout(jday)=dmc
                 dcout(jday)=dc

 290             continue

                  return
                  end
 
         
          
         

       
         
       
       
       
       
       
       
       
       
       
       
     
