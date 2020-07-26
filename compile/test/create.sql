create database stu;

create table info (name CHAR(6), age int);
insert into info values('fff',0);
insert into info values('dec',10);
insert into info values('oxa',16);
insert into info values('hex',16);
insert into info values('uuz',1000);
;;;


create table TeSt (tt iNt);
select * from __dbmeta
select * from __tbmeta_info
sElect * from __tbmeta_test
select * froM info where age<30
select * from info where name='fff' and 'fff'=name

select * from info where age>10
select * from info where age<=10
select * from info where age>10 and age<30
select * frOm INFO where age>10 or age<30
select * from info where name!='fff'
select * from info where name<>'uuz' and name!='fff'

select * from info WHERE age>10 or age<30 and name<>'fff' -- priority AND > OR
select * from info WHERE (age>10 or age<30) and name<>'hex'

create table score (name CHAR(6), score int)
insert into score values('fff',0)
insert into score values('fff',1)
insert into score values('fff',2)
insert into score values('fff',3)
insert into score values('fff',4)
insert into score values('oxa',99)
insert into score values('oxa',100)
insert into score values('ss'); -- field number error
insert into score values(88,100); --value type mismatch

update score set score = 10 where score = 0
update score set score = 0 where score = 10

insert into test (tt) values (4);
insert into test (tt) values (6);
insert into test (tt) values (0);
insert into test (ts) values (6); -- field name error
insert into tast (tt) values (7); -- table name error
insert into test (tt) values ( '6' ) -- value type mismatch

select * from info,score where score<3 and age<5 -- 3 items
select * from info, score -- 35 items
select * from info, score where name='fff' -- name conflict error
select * from info, score where info.name=score.name -- 7 items
select * from info, score where info.name=score.name and (score>2 and score<=99) -- 3 items

delete from test where tt<5;
delete from test
insert into test (tt) values (4);
drop table test;

save


exit

