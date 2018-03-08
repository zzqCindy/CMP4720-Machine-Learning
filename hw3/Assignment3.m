
clear all

global pa pb pxab pcx pdx

pa = [0.25;0.25;0.25;0.25];
pb = [0.6;0.4];
pxab(:,:,1)=[0.5 0.6 0.4 0.2; 0.5 0.4 0.6 0.8];
pxab(:,:,2)=[0.7 0.8 0.1 0.3; 0.3 0.2 0.9 0.7];
pcx=[0.6 0.2; 0.2 0.3; 0.2 0.5];
pdx=[0.3 0.6; 0.7 0.4];

p1 = pabxcd_cond_abxcd('summer','north','sea','dark','thin');
p2 = pabxcd_cond_abxcd('x1','c1','b2');
p3 = pabxcd_cond_abxcd('x2','c1','b2');
p4 = pabxcd_cond_abxcd('south','light');
p5 = pabxcd_cond_abxcd('x1','|','c1','b2');
p6 = pabxcd_cond_abxcd('x2','|','c1','b2');

p2a = pabxcd_cond_abxcd('salmon','|','light','thin','south','winter');
p2b1 = pabxcd_cond_abxcd('winter','|','thin','dark','south');
p2b2 = pabxcd_cond_abxcd('spring','|','thin','dark','south');
p2b3 = pabxcd_cond_abxcd('summer','|','thin','dark','south');
p2b4 = pabxcd_cond_abxcd('autumn','|','thin','dark','south');
p2c = pabxcd_cond_abxcd('north','|','dark','wide','summer');

disp(['P(a3,b1,x2,c3,d2) = ',num2str(p1)]);
disp(['P(x1,c1,b2) = ',num2str(p2)]);
disp(['P(x2,c1,b2) = ',num2str(p3)]);
disp(['P(c1,b2) = ',num2str(p4)]);
disp(['P(x1|c1,b2) = ',num2str(p5)]);
disp(['P(x2|c1,b2) = ',num2str(p6)]);
disp(['P(salmon|light,thin,south,winter) = ',num2str(p2a)]);
disp(['P(winter|thin,dark,south) = ',num2str(p2b1)]);
disp(['P(spring|thin,dark,south) = ',num2str(p2b2)]);
disp(['P(summer|thin,dark,south) = ',num2str(p2b3)]);
disp(['P(autumn|thin,dark,south) = ',num2str(p2b4)]);
disp(['P(north|dark,wide,summer) = ',num2str(p2c)]);

function [a,b,x,c,d,ca,cb,cx,cc,cd] = scaninput(varargin)
    
    flag = 0; a = 0; b = 0; x = 0; c = 0; d = 0;
    ca = 0; cb = 0; cx = 0; cc = 0; cd = 0;
    for v= 1:numel(varargin{1})
    	str = varargin{1}{v};
        switch(str)
            case 'winter'
                s = 'a1';
            case 'spring'
                s = 'a2';
            case 'summer'
                s = 'a3';
            case 'autumn'
                s = 'a4';
            case 'north'
                s = 'b1';
            case 'south'
                s = 'b2';
            case 'salmon'
                s = 'x1';
            case 'sea'
                s = 'x2';
            case 'light'
                s = 'c1';
            case 'medium'
                s = 'c2';
            case 'dark'
                s = 'c3';
            case 'wide'
                s = 'd1';
            case 'thin'
                s = 'd2';
            otherwise
                s = str;
        end
        if flag == 0
            if s(1) == 'a'
                a = str2num(s(2));
            elseif s(1) == 'b'
                b = str2num(s(2));
            elseif s(1) == 'c'
                c = str2num(s(2));
            elseif s(1) == 'd'
                d = str2num(s(2));
            elseif s(1) == 'x'
                x = str2num(s(2));
            else
                flag = 1;
            end 
        end
        if flag == 1
            if s(1) == 'a'
                ca = str2num(s(2));
                a = ca;
            elseif s(1) == 'b'
                cb = str2num(s(2));
                b = cb;
            elseif s(1) == 'c'
                cc = str2num(s(2));
                c = cc;
            elseif s(1) == 'd'
                cd = str2num(s(2));
                d = cd;
            elseif s(1) == 'x'
                cx = str2num(s(2));
                x = cx;
            end
        end
    end

end

function ret = pabxcd(a,b,x,c,d)
    global pa pb pxab pcx pdx
    ret = 1;
    if a ~= 0
        ret = ret * pa(a,1);
    end
    if b ~= 0
        ret = ret * pb(b,1);
    end
    if x ~= 0
        if a ==0 && b == 0
            p = 0;
            for a = 1:4
                for b = 1:2
                    p = p + pxab(x,a,b)*pa(a,1)*pb(b,1);
                end
            end
            ret = ret * p;
        elseif a == 0
            ret = ret * (pxab(x,1,b)*pa(1,1) + pxab(x,2,b)*pa(2,1) + ...
                pxab(x,3,b)*pa(3,1) + pxab(x,4,b)*pa(4,1));
        elseif b == 0
            ret = ret * (pxab(x,a,1)*pb(1,1) + pxab(x,a,2)*pb(2,1));
        else
            ret = ret*pxab(x,a,b);
        end
        if c ~= 0
            ret = ret*pcx(c,x); 
        end
        if d ~= 0
            ret = ret*pdx(d,x);
        end
    else
        mul1 = 1; mul2 = 1;
        if c ~= 0
            mul1 = mul1 * pcx(c,1);
            mul2 = mul2 * pcx(c,2);
        end
        if d ~= 0
            mul1 = mul1 * pdx(d,1);
            mul2 = mul2 * pdx(d,2);
        end
        if c ~= 0 || d ~= 0
            ret = ret*(mul1*pabxcd(a,b,1,0,0) + mul2*pabxcd(a,b,2,0,0));
            if a ~= 0
                ret = ret / pa(a,1);
            end
            if b ~= 0
                ret = ret / pb(b,1);
            end
        end
    end
end

function ret = pabxcd_cond_abxcd(varargin)

    [a,b,x,c,d,ca,cb,cx,cc,cd] = scaninput(varargin);
    p1 = pabxcd(a,b,x,c,d);
    p2 = pabxcd(ca,cb,cx,cc,cd);
    ret = p1/p2;
    
end