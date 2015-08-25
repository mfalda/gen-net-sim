#1 "output\\pars.c"
#1 "<built-in>"
#1 "<command line>"
#1 "output\\pars.c"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<glib.h>
#26 "output\\pars.c"
enum l124{l140,l116,l80,l95};enum l138{l59,l25,l75,l4};enum l126{l53,l19
,l36,l39,l27,l38,l29,l33,l40};enum l149{l96,l90,l83,l68,l76,l99,l88,
l94,l112,l87,l91,l101,l102,l107,l111,l71,l79,l115,l113,l109,l77,l70,
l110,l78,l108,l92,l103,l69,l89,l65,l100,l82,l74};struct lb{GString*ld;double li;
int lp;int l23;enum l138 ls;int lg;struct lb*lz;unsigned l18;};GQueue*l5;GHashTable*l24;
GHashTable*l30;GHashTable*lo;struct lb l17[21+20];struct lb l147[3+7];struct lb l167[2+8];int
l128,l93,l130;int l6[9];void Ripristina();void Cancella();struct lb*l64(const char*ld,int ls);
double l137(const char*ld);double l139(const char*ld);void DefC(const char*ld,double li);void DefV
(const char*ld);GString*l119(const char*ld,int l45);GString*Informa();GString*DefF(const char*ld,const
char*l141,int le,int l45);double l67(int l49);GString*Calcola(const char*ld,double*l105,double*l14);
void Ripristina(){struct lb*l28;struct lb*l12;struct lb*la;l128=33;l93=3;l130=2;
#56 "output\\pars.c"
l6[l53]=(1<<(l53+1))|(1<<(l19+1))|(1<<(l27+1))|(1<<(l29+1))|(1<<(l33+
1));l6[l19]=(1<<(l36+1))|(1<<(l39+1))|(1<<(l38+1))|(1<<(l40+1));l6[
l36]=(1<<(l19+1))|(1<<(l27+1))|(1<<(l29+1))|(1<<(l33+1));l6[l39]=(1<<
(l36+1))|(1<<(l39+1))|(1<<(l38+1))|(1<<(l40+1));l6[l27]=(1<<(l19+1))|
(1<<(l27+1))|(1<<(l29+1))|(1<<(l33+1));l6[l38]=(1<<(l36+1))|(1<<(l39+
1))|(1<<(l38+1))|(1<<(l40+1));l6[l29]=(1<<(l19+1))|(1<<(l27+1))|(1<<(
l29+1))|(1<<(l33+1));l6[l33]=(1<<(l27+1));l6[l40]=(1<<(l19+1))|(1<<(
l27+1))|(1<<(l29+1))|(1<<(l33+1));l30=g_hash_table_new(g_str_hash,g_str_equal);l28=(struct lb* )malloc(1*
sizeof(struct lb));l28->ld=g_string_new("Pi");l28->li=3.14159265358979323846;l28->lg=0;
g_hash_table_insert(l30,"Pi",(gpointer)l28);l28=(struct lb* )malloc(1*sizeof(struct lb));l28->ld=g_string_new("e");l28
->li=exp(1);l28->lg=0;g_hash_table_insert(l30,"e",(gpointer)l28);l24=g_hash_table_new(g_str_hash,g_str_equal);l12=(struct lb
 * )malloc(1*sizeof(struct lb));l12->ld=g_string_new("x");l12->l23=1;l12->lg=0;g_hash_table_insert(l24,"x",(
gpointer)l12);l12=(struct lb* )malloc(1*sizeof(struct lb));l12->ld=g_string_new("y");l12->l23=2;l12->
lg=0;g_hash_table_insert(l24,"y",(gpointer)l12);l12=(struct lb* )malloc(1*sizeof(struct lb));l12->ld=g_string_new("z");
l12->l23=3;l12->lg=0;g_hash_table_insert(l24,"z",(gpointer)l12);lo=g_hash_table_new(g_str_hash,g_str_equal);la=(struct lb* )malloc
(1*sizeof(struct lb));la->ld=g_string_new("?");la->lg=1;la->li=l96;la->ls=l59;la->lp=1;
g_hash_table_insert(lo,"§",(gpointer)la);la=malloc(1*sizeof(struct lb));la->ld=g_string_new("+");la->lg=2;la->li=
l90;la->ls=l25;la->lp=3;g_hash_table_insert(lo,"+",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));
la->ld=g_string_new("-");la->lg=2;la->li=l83;la->ls=l25;la->lp=3;g_hash_table_insert(lo,"-",(gpointer)la
);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("*");la->lg=2;la->li=l68;la->
ls=l25;la->lp=2;g_hash_table_insert(lo,"*",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=
g_string_new("/");la->lg=2;la->li=l76;la->ls=l25;la->lp=2;g_hash_table_insert(lo,"/",(gpointer)la);la=
(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("^");la->lg=2;la->li=l99;la->ls=
l25;la->lp=1;g_hash_table_insert(lo,"^",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new(
"abs");la->lg=1;la->li=l88;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"abs",(gpointer)la);la=
(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("sin");la->lg=1;la->li=l94;la->ls=
l4;la->lp=1;g_hash_table_insert(lo,"sin",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new
("cos");la->lg=1;la->li=l112;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"cos",(gpointer)la);
la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("log");la->lg=2;la->li=l87;la->
ls=l4;la->lp=1;g_hash_table_insert(lo,"log",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld
=g_string_new("atan");la->lg=1;la->li=l91;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"atan",(gpointer)la
);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("tan");la->lg=1;la->li=l101;
la->ls=l4;la->lp=1;g_hash_table_insert(lo,"tan",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la
->ld=g_string_new("%");la->lg=2;la->li=l102;la->ls=l25;la->lp=2;g_hash_table_insert(lo,"%",(gpointer)la
);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("if");la->lg=3;la->li=l107;la
->ls=l4;la->lp=1;g_hash_table_insert(lo,"if",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->
ld=g_string_new("asin");la->lg=1;la->li=l111;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"asin",(
gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("acos");la->lg=1;la->li=
l71;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"acos",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb
));la->ld=g_string_new("!");la->lg=1;la->li=l79;la->ls=l75;la->lp=1;g_hash_table_insert(lo,"!",(
gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("=");la->lg=2;la->li=
l115;la->ls=l25;la->lp=1;g_hash_table_insert(lo,"=",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));
la->ld=g_string_new("~");la->lg=1;la->li=l113;la->ls=l59;la->lp=1;g_hash_table_insert(lo,"~",(gpointer
)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new(">");la->lg=2;la->li=l109;
la->ls=l25;la->lp=1;g_hash_table_insert(lo,">",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la
->ld=g_string_new("<");la->lg=2;la->li=l77;la->ls=l25;la->lp=1;g_hash_table_insert(lo,"<",(gpointer)la
);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("sinh");la->lg=1;la->li=l70;
la->ls=l4;la->lp=1;g_hash_table_insert(lo,"sinh",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));
la->ld=g_string_new("cosh");la->lg=1;la->li=l110;la->ls=l4;la->lp=1;g_hash_table_insert(lo,
"cosh",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("tanh");la->lg=1
;la->li=l78;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"tanh",(gpointer)la);la=(struct lb* )malloc(1*
sizeof(struct lb));la->ld=g_string_new("asinh");la->lg=1;la->li=l108;la->ls=l4;la->lp=1
;g_hash_table_insert(lo,"asinh",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("acosh");
la->lg=1;la->li=l92;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"acosh",(gpointer)la);la=(struct lb
 * )malloc(1*sizeof(struct lb));la->ld=g_string_new("atanh");la->lg=1;la->li=l103;la->ls=l4
;la->lp=1;g_hash_table_insert(lo,"atanh",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new
("sqrt");la->lg=1;la->li=l69;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"sqrt",(gpointer)la);
la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("angle");la->lg=2;la->li=l89;la
->ls=l4;la->lp=1;g_hash_table_insert(lo,"angle",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la
->ld=g_string_new("malloc");la->lg=1;la->li=l65;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"malloc",(gpointer)la
);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("sgn");la->lg=1;la->li=l100;
la->ls=l4;la->lp=1;g_hash_table_insert(lo,"sgn",(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la
->ld=g_string_new("round");la->lg=1;la->li=l82;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"round"
,(gpointer)la);la=(struct lb* )malloc(1*sizeof(struct lb));la->ld=g_string_new("trunc");la->lg=1;la->
li=l74;la->ls=l4;la->lp=1;g_hash_table_insert(lo,"trunc",(gpointer)la);}static void l63(gpointer l168
,gpointer l106,gpointer l156){if(((struct lb* )l106)->ld!=NULL)g_string_free(((struct lb* )l106)->ld,
TRUE);}void Cancella(){g_hash_table_foreach(l30,l63,NULL);g_hash_table_foreach(l24,l63,NULL);g_hash_table_foreach(lo,l63,NULL);}struct lb
 *l64(const char*ld,int ls){struct lb*l14=NULL;if(ls|1)l14=(struct lb* )g_hash_table_lookup(lo,ld);if
(ls|2)l14=(struct lb* )g_hash_table_lookup(l30,ld);if(ls|4)l14=(struct lb* )g_hash_table_lookup(l24,ld);return l14
;}double l137(const char*ld){struct lb*l14=(struct lb* )g_hash_table_lookup(l24,ld);if(l14!=NULL)return l14
->li;else return 0.0;}double l139(const char*ld){struct lb*l14=(struct lb* )g_hash_table_lookup(l30,ld);if(
l14!=NULL)return l14->li;else return 0.0;}void DefC(const char*ld,double li){struct lb*l35;l35
=g_hash_table_lookup(l30,ld);if(l35==NULL){l35=(struct lb* )malloc(1*sizeof(struct lb));l35->ld=g_string_new(ld);
l35->lg=0;g_hash_table_insert(l30,ld,(gpointer)l35);}l35->li=li;return;}void DefV(const char*ld){struct lb
 *l37;l37=g_hash_table_lookup(l24,ld);if(l37==NULL){l37=(struct lb* )malloc(1*sizeof(struct lb));l37->
ld=g_string_new(ld);l37->l23=l93+1;l37->lg=0;g_hash_table_insert(l24,ld,(gpointer)l37);}return;}GString*l119(
const char*ld,int l45){struct lb*la=l64(ld,1);GString*lk=g_string_new("");if(la==NULL){g_string_printf(lk,
"function '%s' does not exists!",ld);return lk;}la->lp=l45;g_string_free(lk,TRUE);return
NULL;}void l131(void*l26,void*l121,GList* *l84){ *l84=g_list_append( *l84,l121);}GString
 *Informa(){GString*l46=g_string_new("");GList*l61=NULL, *l41=NULL;struct lb*l55=NULL;g_string_printf(l46,
"Predefined constants: Pi, e\n");g_string_append_printf(l46,
"Predefined variables: x, y, z\n");g_string_append_printf(l46,
"Predefined functions: name/args (priority)\n");g_hash_table_foreach(lo,l131,&l61);
for(l41=l61;l41!=NULL&&l41->data!=NULL;l41=l41->next){l55=(struct lb* )l41->
data;g_string_append_printf(l46,"\t%s/%d (%d)\n",l55->ld->str,l55->lg,l55->lp);}g_list_free(
l61);return l46;}GString*DefF(const char*ld,const char*ly,int le,int l45){int lf=0,l2=0,
l18;int l60,l21=0;int l54,l43=0;GString*l31=g_string_new("");GString*lk=g_string_new("");struct lb*l3, *
l15, *lx;enum l124 l42=l140;enum l126 l1=l53;lx=g_hash_table_lookup(lo,ld);if(lx==NULL){
lx=(struct lb* )malloc(1*sizeof(struct lb));g_hash_table_insert(lo,ld,(gpointer)lx);lx->lz=NULL;}if(lx->lz!=NULL
){if(lx->lz!=NULL)free(lx->lz);lx->lz=NULL;}l18=strlen(ly);lx->lz=(struct lb* )malloc
(l18*sizeof(struct lb));lx->l18=l18;lx->ld=g_string_new(ld);lx->lg=le;lx->lp=l45;l5=
g_queue_new();while(lf<l18){while(ly[lf]==' ')lf++;if(ly[lf]=='\0')break;else if(g_ascii_isdigit
(ly[lf])){if((l6[l1]&(1<<(l19+1)))==0){g_string_printf(lk,
"syntax error in the expression '%s': unexpected number at position %d!"
,ly,lf+1);Cancella();return lk;}l1=l19;lx->lz[l2].lg=0;lx->lz[l2].ld=g_string_new("");lx
->lz[l2].li=0;l54=0;l60=10;do{if(ly[lf]=='.'){lf++;if(l54==1){g_string_printf(
lk,"syntax error: double decimal at position %d!",lf);Cancella();return lk;}
l54=1;}if(l54==1){lx->lz[l2].li+=(double)(ly[lf]-48)/l60;l60*=10;}else{lx
->lz[l2].li=lx->lz[l2].li*10+(ly[lf]-48);}if(lf>l18)break;lf++;}while(g_ascii_isdigit(
ly[lf])||ly[lf]=='.');l2++;l21++;}else if(ly[lf]=='('){if((l6[l1]&(1<<
(l27+1)))==0){g_string_printf(lk,
"syntax error in the expression '%s': unexpected open parenthesis at position %d!"
,ly,lf+1);Cancella();return lk;}l1=l27;l43++;lf++;}else if(ly[lf]==')'||ly[lf]
==','){if(ly[lf]==')'){if((l6[l1]&(1<<(l38+1)))==0){g_string_printf(lk,
"syntax error in the expression '%s': unexpected closed parenthesis at position %d!"
,ly,lf+1);Cancella();return lk;}l1=l38;}else{if((l6[l1]&l40)==0){g_string_printf(lk,
"syntax error in the expression '%s': unexpected separator at position %d!"
,ly,lf+1);Cancella();return lk;}l1=l40;}l3=(struct lb* )g_queue_peek_tail(l5);while(l3!=NULL&&l3->
l23>=l43&&l3->lg>0){l21-=l3->lg;l21++;lx->lz[l2]= * (struct lb* )g_queue_pop_tail(l5);
l2++;l3=(struct lb* )g_queue_peek_tail(l5);}if(ly[lf]==')')l43--;lf++;}else if(!g_ascii_isdigit(ly[
lf])&&ly[lf]!='('&&ly[lf]!=')'&&ly[lf]!=','){g_string_assign(l31,"");do{g_string_append_c(
l31,ly[lf]);lf++;}while(lf<l18&&g_ascii_isalpha(ly[lf-1])&&g_ascii_isalpha(ly[lf]));if((l1==
l53||l1==l29||l1==l27||l1==l36)&&ly[lf-1]=='-')g_string_assign(l31,"§");if((l15=
g_hash_table_lookup(lo,l31->str))!=NULL)l42=l95;else if((l15=g_hash_table_lookup(l30,l31->str))!=NULL)l42=
l116;else if((l15=g_hash_table_lookup(l24,l31->str))!=NULL)l42=l80;if(l15==NULL){g_string_printf(lk,
"syntax error in the expression '%s': symbol '%s' does not exists",ly
,l31->str);Cancella();return lk;}if(l42==l116){if((l6[l1]&(1<<(l19+1)))==0){
g_string_printf(lk,
"syntax error in the expression '%s': unexpected constant at position %d!"
,ly,lf);Cancella();return lk;}l1=l19;lx->lz[l2].ld=g_string_new(l15->ld->str);lx->lz[l2]
.lg=0;lx->lz[l2].li=l15->li;l2++;l21++;}else if(l42==l80){if((l6[l1]&(
1<<(l19+1)))==0){g_string_printf(lk,
"syntax error in the expression '%s': unexpected variable at position %d!"
,ly,lf);Cancella();return lk;}l1=l19;lx->lz[l2].ld=NULL;lx->lz[l2].lg=0;lx->lz[
l2].l23=l15->l23;l2++;l21++;}else if(l42==l95){switch(l15->ls){case l59:if
((l6[l1]&(1<<(l29+1)))==0){g_string_printf(lk,
"syntax error: unexpected prefix operator at position %d!",lf);Cancella();
return lk;}l1=l29;break;case l4:if((l6[l1]&(1<<(l33+1)))==0){g_string_printf(lk,
"syntax error: unexpected function at position %d!",lf);Cancella();return lk;}
l1=l33;break;case l25:if((l6[l1]&(1<<(l36+1)))==0){g_string_printf(lk,
"syntax error: unexpected infix operator at position %d!",lf);Cancella();
return lk;}l1=l36;break;case l75:if((l6[l1]&(1<<(l39+1)))==0){g_string_printf(lk,
"syntax error: unexpected postfix operator at position %d!",lf);Cancella();
return lk;}l1=l39;break;}l3=(struct lb* )g_queue_peek_tail(l5);l15->l23=l43;if(l15!=NULL&&(l3==
NULL||(l15->lp>=l3->lp&&l43<=l3->l23))){while(l3!=NULL&&l3->lp<=l15->lp){l3
=(struct lb* )g_queue_pop_tail(l5);l21-=l3->lg;lx->lz[l2]= *l3;l2++;l21++;l3=(struct lb* )g_queue_peek_tail
(l5);}}g_queue_push_tail(l5,l15);}}}l3=(struct lb* )g_queue_pop_tail(l5);while(l3!=NULL&&l3->lp>0){l21
-=l3->lg;lx->lz[l2]= *l3;l2++;l21++;l3=(struct lb* )g_queue_pop_tail(l5);}if(l21<1){
g_string_assign(lk,"missing arguments!");Cancella();return lk;}if(l21>1){g_string_assign(lk,
"too much arguments!");Cancella();return lk;}if(l43!=0){g_string_assign(lk,
"parentheses do not match!");Cancella();return lk;}lx->l18=l2;g_queue_free(l5);g_string_free(l31
,TRUE);g_string_free(lk,TRUE);return NULL;}double l67(int l49){int l26;if(l49==0)return 1;for(
l26=l49-1;l26>0;l26--)l49*=l26;return l49;}GString*Calcola(const char*ld,double*l105,double*
l14){double le[10];GString*l62=g_string_new("");GString*lk=g_string_new("");int l85,l51=0;double lv, *l5;
int l26;unsigned int lf=0;struct lb*l17;l17=(struct lb* )g_hash_table_lookup(lo,ld);if(l17==NULL){g_string_printf
(lk,"semantic error: function '%s' does not exists!",ld);return lk;}l5=(
double* )malloc(l17->l18*sizeof(double));while(lf<l17->l18){if(l17->lz[lf].lg==0){if(
l17->lz[lf].ld!=NULL){l5[l51++]=l17->lz[lf].li;}else{l5[l51++]=l105[l17
->lz[lf].l23-1];}lf++;}else{for(l26=l17->lz[lf].lg-1;l26>-1;l26--)le[
l26]=l5[--l51];l85=(int)l17->lz[lf].li;g_string_assign(l62,l17->lz[lf].ld->str);
switch(l85){case l96:lv=-le[0];break;case l90:lv=le[0]+le[1];break;case l83:lv=le[0
]-le[1];break;case l68:lv=le[0] *le[1];break;case l76:lv=le[0]/le[1];break;case l99:
lv=pow(le[0],le[1]);break;case l88:lv=fabs(le[0]);break;case l94:lv=sin(le[0]
);break;case l112:lv=cos(le[0]);break;case l87:lv=log(le[0])/log(le[1]);break;case
l91:lv=atan(le[0]);break;case l101:lv=tan(le[0]);break;case l102:lv=fmod(le[0]
,le[1]);break;case l107:lv=(le[0]>0.0)?le[1]:le[2];break;case l111:lv=asin(le[0
]);break;case l71:lv=acos(le[0]);break;case l79:lv=l67((int)le[0]);break;case l115:lv
=fabs(le[0]-le[1])<=DBL_EPSILON*fabs(le[0]);break;case l113:lv=!(le[0]>0.0);break;case
l77:lv=le[0]>le[1];break;case l109:lv=le[0]<le[1];break;case l70:lv=sinh(le[0]);
break;case l110:lv=cosh(le[0]);break;case l78:lv=tanh(le[0]);break;case l108:lv=log(
le[0]+sqrt(le[0] *le[0]+1));break;case l92:lv=log(le[0]+sqrt(le[0] *le[0]-1));
break;case l103:lv=log((1+le[0])/(1-le[0]))/2;break;case l69:lv=sqrt(le[0]);break;
case l89:lv=atan2(le[0],le[1]);break;case l65:lv=log(le[0]);break;case l100:lv=(le
[0]>=0)?((le[0]>0)?1:0):-1;break;case l82:lv=floor(le[0]+0.5);break;case l74:lv=
(int)le[0];break;default:lk=Calcola(l62->str,le,&lv);if(lk!=NULL){g_string_printf(lk,
"error: %s!",lk->str);return lk;}}l5[l51++]=lv;if(lf==l17->l18)break;lf++;}}
 *l14=l5[--l51];if(l5!=NULL)free(l5);g_string_free(l62,TRUE);g_string_free(lk,TRUE);return NULL;}
