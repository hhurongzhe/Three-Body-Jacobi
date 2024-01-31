// MIT License

// Copyright (c) 2022 0382

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once
#ifndef JSHL_WIGNERSYMBOL_HPP
#define JSHL_WIGNERSYMBOL_HPP

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace util
{

    // clang-format off
constexpr uint64_t u64_binomial_data[] = {1ull,1ull,1ull,2ull,1ull,3ull,1ull,4ull,6ull,1ull,5ull,10ull,1ull,6ull,15ull,20ull,1ull,7ull,21ull,35ull,1ull,8ull,28ull,56ull,70ull,1ull,9ull,36ull,84ull,126ull,1ull,10ull,45ull,120ull,210ull,252ull,1ull,11ull,55ull,165ull,330ull,462ull,1ull,12ull,66ull,220ull,495ull,792ull,924ull,1ull,13ull,78ull,286ull,715ull,1287ull,1716ull,1ull,14ull,91ull,364ull,1001ull,2002ull,3003ull,3432ull,1ull,15ull,105ull,455ull,1365ull,3003ull,5005ull,6435ull,1ull,16ull,120ull,560ull,1820ull,4368ull,8008ull,11440ull,12870ull,1ull,17ull,136ull,680ull,2380ull,6188ull,12376ull,19448ull,24310ull,1ull,18ull,153ull,816ull,3060ull,8568ull,18564ull,31824ull,43758ull,48620ull,1ull,19ull,171ull,969ull,3876ull,11628ull,27132ull,50388ull,75582ull,92378ull,1ull,20ull,190ull,1140ull,4845ull,15504ull,38760ull,77520ull,125970ull,167960ull,184756ull,1ull,21ull,210ull,1330ull,5985ull,20349ull,54264ull,116280ull,203490ull,293930ull,352716ull,1ull,22ull,231ull,1540ull,7315ull,26334ull,74613ull,170544ull,319770ull,497420ull,646646ull,705432ull,1ull,23ull,253ull,1771ull,8855ull,33649ull,100947ull,245157ull,490314ull,817190ull,1144066ull,1352078ull,1ull,24ull,276ull,2024ull,10626ull,42504ull,134596ull,346104ull,735471ull,1307504ull,1961256ull,2496144ull,2704156ull,1ull,25ull,300ull,2300ull,12650ull,53130ull,177100ull,480700ull,1081575ull,2042975ull,3268760ull,4457400ull,5200300ull,1ull,26ull,325ull,2600ull,14950ull,65780ull,230230ull,657800ull,1562275ull,3124550ull,5311735ull,7726160ull,9657700ull,10400600ull,1ull,27ull,351ull,2925ull,17550ull,80730ull,296010ull,888030ull,2220075ull,4686825ull,8436285ull,13037895ull,17383860ull,20058300ull,1ull,28ull,378ull,3276ull,20475ull,98280ull,376740ull,1184040ull,3108105ull,6906900ull,13123110ull,21474180ull,30421755ull,37442160ull,40116600ull,1ull,29ull,406ull,3654ull,23751ull,118755ull,475020ull,1560780ull,4292145ull,10015005ull,20030010ull,34597290ull,51895935ull,67863915ull,77558760ull,1ull,30ull,435ull,4060ull,27405ull,142506ull,593775ull,2035800ull,5852925ull,14307150ull,30045015ull,54627300ull,86493225ull,119759850ull,145422675ull,155117520ull,1ull,31ull,465ull,4495ull,31465ull,169911ull,736281ull,2629575ull,7888725ull,20160075ull,44352165ull,84672315ull,141120525ull,206253075ull,265182525ull,300540195ull,1ull,32ull,496ull,4960ull,35960ull,201376ull,906192ull,3365856ull,10518300ull,28048800ull,64512240ull,129024480ull,225792840ull,347373600ull,471435600ull,565722720ull,601080390ull,1ull,33ull,528ull,5456ull,40920ull,237336ull,1107568ull,4272048ull,13884156ull,38567100ull,92561040ull,193536720ull,354817320ull,573166440ull,818809200ull,1037158320ull,1166803110ull,1ull,34ull,561ull,5984ull,46376ull,278256ull,1344904ull,5379616ull,18156204ull,52451256ull,131128140ull,286097760ull,548354040ull,927983760ull,1391975640ull,1855967520ull,2203961430ull,2333606220ull,1ull,35ull,595ull,6545ull,52360ull,324632ull,1623160ull,6724520ull,23535820ull,70607460ull,183579396ull,417225900ull,834451800ull,1476337800ull,2319959400ull,3247943160ull,4059928950ull,4537567650ull,1ull,36ull,630ull,7140ull,58905ull,376992ull,1947792ull,8347680ull,30260340ull,94143280ull,254186856ull,600805296ull,1251677700ull,2310789600ull,3796297200ull,5567902560ull,7307872110ull,8597496600ull,9075135300ull,1ull,37ull,666ull,7770ull,66045ull,435897ull,2324784ull,10295472ull,38608020ull,124403620ull,348330136ull,854992152ull,1852482996ull,3562467300ull,6107086800ull,9364199760ull,12875774670ull,15905368710ull,17672631900ull,1ull,38ull,703ull,8436ull,73815ull,501942ull,2760681ull,12620256ull,48903492ull,163011640ull,472733756ull,1203322288ull,2707475148ull,5414950296ull,9669554100ull,15471286560ull,22239974430ull,28781143380ull,33578000610ull,35345263800ull,1ull,39ull,741ull,9139ull,82251ull,575757ull,3262623ull,15380937ull,61523748ull,211915132ull,635745396ull,1676056044ull,3910797436ull,8122425444ull,15084504396ull,25140840660ull,37711260990ull,51021117810ull,62359143990ull,68923264410ull,1ull,40ull,780ull,9880ull,91390ull,658008ull,3838380ull,18643560ull,76904685ull,273438880ull,847660528ull,2311801440ull,5586853480ull,12033222880ull,23206929840ull,40225345056ull,62852101650ull,88732378800ull,113380261800ull,131282408400ull,137846528820ull,1ull,41ull,820ull,10660ull,101270ull,749398ull,4496388ull,22481940ull,95548245ull,350343565ull,1121099408ull,3159461968ull,7898654920ull,17620076360ull,35240152720ull,63432274896ull,103077446706ull,151584480450ull,202112640600ull,244662670200ull,269128937220ull,1ull,42ull,861ull,11480ull,111930ull,850668ull,5245786ull,26978328ull,118030185ull,445891810ull,1471442973ull,4280561376ull,11058116888ull,25518731280ull,52860229080ull,98672427616ull,166509721602ull,254661927156ull,353697121050ull,446775310800ull,513791607420ull,538257874440ull,1ull,43ull,903ull,12341ull,123410ull,962598ull,6096454ull,32224114ull,145008513ull,563921995ull,1917334783ull,5752004349ull,15338678264ull,36576848168ull,78378960360ull,151532656696ull,265182149218ull,421171648758ull,608359048206ull,800472431850ull,960566918220ull,1052049481860ull,1ull,44ull,946ull,13244ull,135751ull,1086008ull,7059052ull,38320568ull,177232627ull,708930508ull,2481256778ull,7669339132ull,21090682613ull,51915526432ull,114955808528ull,229911617056ull,416714805914ull,686353797976ull,1029530696964ull,1408831480056ull,1761039350070ull,2012616400080ull,2104098963720ull,1ull,45ull,990ull,14190ull,148995ull,1221759ull,8145060ull,45379620ull,215553195ull,886163135ull,3190187286ull,10150595910ull,28760021745ull,73006209045ull,166871334960ull,344867425584ull,646626422970ull,1103068603890ull,1715884494940ull,2438362177020ull,3169870830126ull,3773655750150ull,4116715363800ull,1ull,46ull,1035ull,15180ull,163185ull,1370754ull,9366819ull,53524680ull,260932815ull,1101716330ull,4076350421ull,13340783196ull,38910617655ull,101766230790ull,239877544005ull,511738760544ull,991493848554ull,1749695026860ull,2818953098830ull,4154246671960ull,5608233007146ull,6943526580276ull,7890371113950ull,8233430727600ull,1ull,47ull,1081ull,16215ull,178365ull,1533939ull,10737573ull,62891499ull,314457495ull,1362649145ull,5178066751ull,17417133617ull,52251400851ull,140676848445ull,341643774795ull,751616304549ull,1503232609098ull,2741188875414ull,4568648125690ull,6973199770790ull,9762479679106ull,12551759587422ull,14833897694226ull,16123801841550ull,1ull,48ull,1128ull,17296ull,194580ull,1712304ull,12271512ull,73629072ull,377348994ull,1677106640ull,6540715896ull,22595200368ull,69668534468ull,192928249296ull,482320623240ull,1093260079344ull,2254848913647ull,4244421484512ull,7309837001104ull,11541847896480ull,16735679449896ull,22314239266528ull,27385657281648ull,30957699535776ull,32247603683100ull,1ull,49ull,1176ull,18424ull,211876ull,1906884ull,13983816ull,85900584ull,450978066ull,2054455634ull,8217822536ull,29135916264ull,92263734836ull,262596783764ull,675248872536ull,1575580702584ull,3348108992991ull,6499270398159ull,11554258485616ull,18851684897584ull,28277527346376ull,39049918716424ull,49699896548176ull,58343356817424ull,63205303218876ull,1ull,50ull,1225ull,19600ull,230300ull,2118760ull,15890700ull,99884400ull,536878650ull,2505433700ull,10272278170ull,37353738800ull,121399651100ull,354860518600ull,937845656300ull,2250829575120ull,4923689695575ull,9847379391150ull,18053528883775ull,30405943383200ull,47129212243960ull,67327446062800ull,88749815264600ull,108043253365600ull,121548660036300ull,126410606437752ull,1ull,51ull,1275ull,20825ull,249900ull,2349060ull,18009460ull,115775100ull,636763050ull,3042312350ull,12777711870ull,47626016970ull,158753389900ull,476260169700ull,1292706174900ull,3188675231420ull,7174519270695ull,14771069086725ull,27900908274925ull,48459472266975ull,77535155627160ull,114456658306760ull,156077261327400ull,196793068630200ull,229591913401900ull,247959266474052ull,1ull,52ull,1326ull,22100ull,270725ull,2598960ull,20358520ull,133784560ull,752538150ull,3679075400ull,15820024220ull,60403728840ull,206379406870ull,635013559600ull,1768966344600ull,4481381406320ull,10363194502115ull,21945588357420ull,42671977361650ull,76360380541900ull,125994627894135ull,191991813933920ull,270533919634160ull,352870329957600ull,426384982032100ull,477551179875952ull,495918532948104ull,1ull,53ull,1378ull,23426ull,292825ull,2869685ull,22957480ull,154143080ull,886322710ull,4431613550ull,19499099620ull,76223753060ull,266783135710ull,841392966470ull,2403979904200ull,6250347750920ull,14844575908435ull,32308782859535ull,64617565719070ull,119032357903550ull,202355008436035ull,317986441828055ull,462525733568080ull,623404249591760ull,779255311989700ull,903936161908052ull,973469712824056ull,1ull,54ull,1431ull,24804ull,316251ull,3162510ull,25827165ull,177100560ull,1040465790ull,5317936260ull,23930713170ull,95722852680ull,343006888770ull,1108176102180ull,3245372870670ull,8654327655120ull,21094923659355ull,47153358767970ull,96926348578605ull,183649923622620ull,321387366339585ull,520341450264090ull,780512175396135ull,1085929983159840ull,1402659561581460ull,1683191473897752ull,1877405874732108ull,1946939425648112ull,1ull,55ull,1485ull,26235ull,341055ull,3478761ull,28989675ull,202927725ull,1217566350ull,6358402050ull,29248649430ull,119653565850ull,438729741450ull,1451182990950ull,4353548972850ull,11899700525790ull,29749251314475ull,68248282427325ull,144079707346575ull,280576272201225ull,505037289962205ull,841728816603675ull,1300853625660225ull,1866442158555975ull,2488589544741300ull,3085851035479212ull,3560597348629860ull,3824345300380220ull,1ull,56ull,1540ull,27720ull,367290ull,3819816ull,32468436ull,231917400ull,1420494075ull,7575968400ull,35607051480ull,148902215280ull,558383307300ull,1889912732400ull,5804731963800ull,16253249498640ull,41648951840265ull,97997533741800ull,212327989773900ull,424655979547800ull,785613562163430ull,1346766106565880ull,2142582442263900ull,3167295784216200ull,4355031703297275ull,5574440580220512ull,6646448384109072ull,7384942649010080ull,7648690600760440ull,1ull,57ull,1596ull,29260ull,395010ull,4187106ull,36288252ull,264385836ull,1652411475ull,8996462475ull,43183019880ull,184509266760ull,707285522580ull,2448296039700ull,7694644696200ull,22057981462440ull,57902201338905ull,139646485582065ull,310325523515700ull,636983969321700ull,1210269541711230ull,2132379668729310ull,3489348548829780ull,5309878226480100ull,7522327487513475ull,9929472283517787ull,12220888964329584ull,14031391033119152ull,15033633249770520ull,1ull,58ull,1653ull,30856ull,424270ull,4582116ull,40475358ull,300674088ull,1916797311ull,10648873950ull,52179482355ull,227692286640ull,891794789340ull,3155581562280ull,10142940735900ull,29752626158640ull,79960182801345ull,197548686920970ull,449972009097765ull,947309492837400ull,1847253511032930ull,3342649210440540ull,5621728217559090ull,8799226775309880ull,12832205713993575ull,17451799771031262ull,22150361247847371ull,26252279997448736ull,29065024282889672ull,30067266499541040ull,1ull,59ull,1711ull,32509ull,455126ull,5006386ull,45057474ull,341149446ull,2217471399ull,12565671261ull,62828356305ull,279871768995ull,1119487075980ull,4047376351620ull,13298522298180ull,39895566894540ull,109712808959985ull,277508869722315ull,647520696018735ull,1397281501935165ull,2794563003870330ull,5189902721473470ull,8964377427999630ull,14420954992868970ull,21631432489303455ull,30284005485024837ull,39602161018878633ull,48402641245296107ull,55317304280338408ull,59132290782430712ull,1ull,60ull,1770ull,34220ull,487635ull,5461512ull,50063860ull,386206920ull,2558620845ull,14783142660ull,75394027566ull,342700125300ull,1399358844975ull,5166863427600ull,17345898649800ull,53194089192720ull,149608375854525ull,387221678682300ull,925029565741050ull,2044802197953900ull,4191844505805495ull,7984465725343800ull,14154280149473100ull,23385332420868600ull,36052387482172425ull,51915437974328292ull,69886166503903470ull,88004802264174740ull,103719945525634515ull,114449595062769120ull,118264581564861424ull,1ull,61ull,1830ull,35990ull,521855ull,5949147ull,55525372ull,436270780ull,2944827765ull,17341763505ull,90177170226ull,418094152866ull,1742058970275ull,6566222272575ull,22512762077400ull,70539987842520ull,202802465047245ull,536830054536825ull,1312251244423350ull,2969831763694950ull,6236646703759395ull,12176310231149295ull,22138745874816900ull,37539612570341700ull,59437719903041025ull,87967825456500717ull,121801604478231762ull,157890968768078210ull,191724747789809255ull,218169540588403635ull,232714176627630544ull,1ull,62ull,1891ull,37820ull,557845ull,6471002ull,61474519ull,491796152ull,3381098545ull,20286591270ull,107518933731ull,508271323092ull,2160153123141ull,8308281242850ull,29078984349975ull,93052749919920ull,273342452889765ull,739632519584070ull,1849081298960175ull,4282083008118300ull,9206478467454345ull,18412956934908690ull,34315056105966195ull,59678358445158600ull,96977332473382725ull,147405545359541742ull,209769429934732479ull,279692573246309972ull,349615716557887465ull,409894288378212890ull,450883717216034179ull,465428353255261088ull,1ull,63ull,1953ull,39711ull,595665ull,7028847ull,67945521ull,553270671ull,3872894697ull,23667689815ull,127805525001ull,615790256823ull,2668424446233ull,10468434365991ull,37387265592825ull,122131734269895ull,366395202809685ull,1012974972473835ull,2588713818544245ull,6131164307078475ull,13488561475572645ull,27619435402363035ull,52728013040874885ull,93993414551124795ull,156655690918541325ull,244382877832924467ull,357174975294274221ull,489462003181042451ull,629308289804197437ull,759510004936100355ull,860778005594247069ull,916312070471295267ull,1ull,64ull,2016ull,41664ull,635376ull,7624512ull,74974368ull,621216192ull,4426165368ull,27540584512ull,151473214816ull,743595781824ull,3284214703056ull,13136858812224ull,47855699958816ull,159518999862720ull,488526937079580ull,1379370175283520ull,3601688791018080ull,8719878125622720ull,19619725782651120ull,41107996877935680ull,80347448443237920ull,146721427591999680ull,250649105469666120ull,401038568751465792ull,601557853127198688ull,846636978475316672ull,1118770292985239888ull,1388818294740297792ull,1620288010530347424ull,1777090076065542336ull,1832624140942590534ull,1ull,65ull,2080ull,43680ull,677040ull,8259888ull,82598880ull,696190560ull,5047381560ull,31966749880ull,179013799328ull,895068996640ull,4027810484880ull,16421073515280ull,60992558771040ull,207374699821536ull,648045936942300ull,1867897112363100ull,4981058966301600ull,12321566916640800ull,28339603908273840ull,60727722660586800ull,121455445321173600ull,227068876035237600ull,397370533061665800ull,651687674221131912ull,1002596421878664480ull,1448194831602515360ull,1965407271460556560ull,2507588587725537680ull,3009106305270645216ull,3397378086595889760ull,3609714217008132870ull,1ull,66ull,2145ull,45760ull,720720ull,8936928ull,90858768ull,778789440ull,5743572120ull,37014131440ull,210980549208ull,1074082795968ull,4922879481520ull,20448884000160ull,77413632286320ull,268367258592576ull,855420636763836ull,2515943049305400ull,6848956078664700ull,17302625882942400ull,40661170824914640ull,89067326568860640ull,182183167981760400ull,348524321356411200ull,624439409096903400ull,1049058207282797712ull,1654284096099796392ull,2450791253481179840ull,3413602103063071920ull,4472995859186094240ull,5516694892996182896ull,6406484391866534976ull,7007092303604022630ull,7219428434016265740ull,1ull,67ull,2211ull,47905ull,766480ull,9657648ull,99795696ull,869648208ull,6522361560ull,42757703560ull,247994680648ull,1285063345176ull,5996962277488ull,25371763481680ull,97862516286480ull,345780890878896ull,1123787895356412ull,3371363686069236ull,9364899127970100ull,24151581961607100ull,57963796707857040ull,129728497393775280ull,271250494550621040ull,530707489338171600ull,972963730453314600ull,1673497616379701112ull,2703342303382594104ull,4105075349580976232ull,5864393356544251760ull,7886597962249166160ull,9989690752182277136ull,11923179284862717872ull,13413576695470557606ull,14226520737620288370ull};
    // clang-format on

    class WignerSymbols
    {
    public:
        WignerSymbols() : _nmax(67)
        {
            // initialize the data
            _binomial_data.resize(_binomial_data_size(_nmax));
            std::copy(std::begin(u64_binomial_data), std::end(u64_binomial_data), _binomial_data.begin());
        }

        // judge if a number is a odd number
        static bool isodd(int x) { return x % 2 != 0; }
        // judge if a number is a even number
        static bool iseven(int x) { return x % 2 == 0; }
        // judge if two number are same odd or same even
        static bool is_same_parity(int x, int y) { return iseven(x ^ y); }
        // return (-1)^n
        static int iphase(int x) { return iseven(x) - isodd(x); }
        // check if m-quantum number if one of the components of a the j-quantum number
        static bool check_jm(int dj, int dm) { return is_same_parity(dj, dm) && (std::abs(dm) <= dj); }
        // judge if three angular momentum can couple
        static bool check_couple(int dj1, int dj2, int dj3)
        {
            return dj1 >= 0 && dj2 >= 0 && is_same_parity(dj1 + dj2, dj3) && (dj3 <= (dj1 + dj2)) &&
                   (dj3 >= std::abs(dj1 - dj2));
        }
        // only works for positive n
        static double quick_pow(double x, int n)
        {
            double ans = 1;
            while (n)
            {
                if (n & 1)
                    ans = ans * x;
                n = n >> 1;
                x = x * x;
            }
            return ans;
        }

        double binomial(int n, int k) const
        {
            if (unsigned(n) > unsigned(_nmax) || unsigned(k) > unsigned(n))
                return 0;
            else
            {
                k = std::min(k, n - k);
                return _binomial_data[_binomial_index(n, k)];
            }
        }

        double unsafe_binomial(int n, int k) const
        {
            k = std::min(k, n - k);
            return _binomial_data[_binomial_index(n, k)];
        }

        double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3) const
        {
            if (!(check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3)))
                return 0;
            if (!check_couple(dj1, dj2, dj3))
                return 0;
            if (dm1 + dm2 != dm3)
                return 0;
            const int J = (dj1 + dj2 + dj3) / 2;
            const int jm1 = J - dj1;
            const int jm2 = J - dj2;
            const int jm3 = J - dj3;
            const int j1mm1 = (dj1 - dm1) / 2;
            const int j2mm2 = (dj2 - dm2) / 2;
            const int j3mm3 = (dj3 - dm3) / 2;
            const int j2pm2 = (dj2 + dm2) / 2;
            const double A = std::sqrt(unsafe_binomial(dj1, jm2) * unsafe_binomial(dj2, jm3) /
                                       (unsafe_binomial(J + 1, jm3) * unsafe_binomial(dj1, j1mm1) *
                                        unsafe_binomial(dj2, j2mm2) * unsafe_binomial(dj3, j3mm3)));
            double B = 0;
            const int low = std::max(0, std::max(j1mm1 - jm2, j2pm2 - jm1));
            const int high = std::min(jm3, std::min(j1mm1, j2pm2));
            for (auto z = low; z <= high; ++z)
            {
                B = -B + unsafe_binomial(jm3, z) * unsafe_binomial(jm2, j1mm1 - z) * unsafe_binomial(jm1, j2pm2 - z);
            }
            return iphase(high) * A * B;
        }

        double CG0(int j1, int j2, int j3) const
        {
            if (!check_couple(2 * j1, 2 * j2, 2 * j3))
                return 0;
            const int J = j1 + j2 + j3;
            if (isodd(J))
                return 0;
            const int g = J / 2;
            return iphase(g - j3) * unsafe_binomial(g, j3) * unsafe_binomial(j3, g - j1) /
                   std::sqrt(unsafe_binomial(J + 1, 2 * j3 + 1) * unsafe_binomial(2 * j3, J - 2 * j1));
        }

        double f3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3) const
        {
            if (!(check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3)))
                return 0;
            if (!check_couple(dj1, dj2, dj3))
                return 0;
            if (dm1 + dm2 + dm3 != 0)
                return 0;
            const int J = (dj1 + dj2 + dj3) / 2;
            const int jm1 = J - dj1;
            const int jm2 = J - dj2;
            const int jm3 = J - dj3;
            const int j1mm1 = (dj1 - dm1) / 2;
            const int j2mm2 = (dj2 - dm2) / 2;
            const int j3mm3 = (dj3 - dm3) / 2;
            const int j1pm1 = (dj1 + dm1) / 2;
            const double A = std::sqrt(unsafe_binomial(dj1, jm2) * unsafe_binomial(dj2, jm1) /
                                       ((J + 1) * unsafe_binomial(J, jm3) * unsafe_binomial(dj1, j1mm1) *
                                        unsafe_binomial(dj2, j2mm2) * unsafe_binomial(dj3, j3mm3)));
            double B = 0;
            const int low = std::max(0, std::max(j1pm1 - jm2, j2mm2 - jm1));
            const int high = std::min(jm3, std::min(j1pm1, j2mm2));
            for (auto z = low; z <= high; ++z)
            {
                B = -B + unsafe_binomial(jm3, z) * unsafe_binomial(jm2, j1pm1 - z) * unsafe_binomial(jm1, j2mm2 - z);
            }
            return iphase(dj1 + (dj3 + dm3) / 2 + high) * A * B;
        }

        double f6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6) const
        {
            if (!(check_couple(dj1, dj2, dj3) && check_couple(dj1, dj5, dj6) && check_couple(dj4, dj2, dj6) &&
                  check_couple(dj4, dj5, dj3)))
                return 0;
            const int j123 = (dj1 + dj2 + dj3) / 2;
            const int j156 = (dj1 + dj5 + dj6) / 2;
            const int j426 = (dj4 + dj2 + dj6) / 2;
            const int j453 = (dj4 + dj5 + dj3) / 2;
            const int jpm123 = (dj1 + dj2 - dj3) / 2;
            const int jpm132 = (dj1 + dj3 - dj2) / 2;
            const int jpm231 = (dj2 + dj3 - dj1) / 2;
            const int jpm156 = (dj1 + dj5 - dj6) / 2;
            const int jpm426 = (dj4 + dj2 - dj6) / 2;
            const int jpm453 = (dj4 + dj5 - dj3) / 2;
            const double A = std::sqrt(unsafe_binomial(j123 + 1, dj1 + 1) * unsafe_binomial(dj1, jpm123) /
                                       (unsafe_binomial(j156 + 1, dj1 + 1) * unsafe_binomial(dj1, jpm156) *
                                        unsafe_binomial(j453 + 1, dj4 + 1) * unsafe_binomial(dj4, jpm453) *
                                        unsafe_binomial(j426 + 1, dj4 + 1) * unsafe_binomial(dj4, jpm426)));
            double B = 0;
            const int low = std::max(j123, std::max(j156, std::max(j426, j453)));
            const int high = std::min(jpm123 + j453, std::min(jpm132 + j426, jpm231 + j156));
            for (auto x = low; x <= high; ++x)
            {
                B = -B + unsafe_binomial(x + 1, j123 + 1) * unsafe_binomial(jpm123, x - j453) *
                             unsafe_binomial(jpm132, x - j426) * unsafe_binomial(jpm231, x - j156);
            }
            return iphase(high) * A * B / (dj4 + 1);
        }

        double Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6) const
        {
            return iphase((dj1 + dj2 + dj3 + dj4) / 2) * f6j(dj1, dj2, dj5, dj4, dj3, dj6);
        }

        double f9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9) const
        {
            if (!(check_couple(dj1, dj2, dj3) && check_couple(dj4, dj5, dj6) && check_couple(dj7, dj8, dj9) &&
                  check_couple(dj1, dj4, dj7) && check_couple(dj2, dj5, dj8) && check_couple(dj3, dj6, dj9)))
                return 0;
            const int j123 = (dj1 + dj2 + dj3) / 2;
            const int j456 = (dj4 + dj5 + dj6) / 2;
            const int j789 = (dj7 + dj8 + dj9) / 2;
            const int j147 = (dj1 + dj4 + dj7) / 2;
            const int j258 = (dj2 + dj5 + dj8) / 2;
            const int j369 = (dj3 + dj6 + dj9) / 2;
            const int pm123 = (dj1 + dj2 - dj3) / 2;
            const int pm132 = (dj1 + dj3 - dj2) / 2;
            const int pm231 = (dj2 + dj3 - dj1) / 2;
            const int pm456 = (dj4 + dj5 - dj6) / 2;
            const int pm465 = (dj4 + dj6 - dj5) / 2;
            const int pm564 = (dj5 + dj6 - dj4) / 2;
            const int pm789 = (dj7 + dj8 - dj9) / 2;
            const int pm798 = (dj7 + dj9 - dj8) / 2;
            const int pm897 = (dj8 + dj9 - dj7) / 2;
            const double P0_nu = unsafe_binomial(j123 + 1, dj1 + 1) * unsafe_binomial(dj1, pm123) * //
                                 unsafe_binomial(j456 + 1, dj5 + 1) * unsafe_binomial(dj5, pm456) * //
                                 unsafe_binomial(j789 + 1, dj9 + 1) * unsafe_binomial(dj9, pm798);
            const double P0_de = unsafe_binomial(j147 + 1, dj1 + 1) * unsafe_binomial(dj1, (dj1 + dj4 - dj7) / 2) *
                                 unsafe_binomial(j258 + 1, dj5 + 1) * unsafe_binomial(dj5, (dj2 + dj5 - dj8) / 2) *
                                 unsafe_binomial(j369 + 1, dj9 + 1) * unsafe_binomial(dj9, (dj3 + dj9 - dj6) / 2);
            const double P0 = std::sqrt(P0_nu / P0_de);
            const int dtl = std::max(std::abs(dj2 - dj6), std::max(std::abs(dj4 - dj8), std::abs(dj1 - dj9)));
            const int dth = std::min(dj2 + dj6, std::min(dj4 + dj8, dj1 + dj9));
            double PABC = 0;
            for (auto dt = dtl; dt <= dth; dt += 2)
            {
                const int j19t = (dj1 + dj9 + dt) / 2;
                const int j26t = (dj2 + dj6 + dt) / 2;
                const int j48t = (dj4 + dj8 + dt) / 2;
                double Pt_de = unsafe_binomial(j19t + 1, dt + 1) * unsafe_binomial(dt, (dj1 + dt - dj9) / 2) *
                               unsafe_binomial(j26t + 1, dt + 1) * unsafe_binomial(dt, (dj2 + dt - dj6) / 2) *
                               unsafe_binomial(j48t + 1, dt + 1) * unsafe_binomial(dt, (dj4 + dt - dj8) / 2);
                Pt_de *= (dt + 1) * (dt + 1);
                const int xl = std::max(j123, std::max(j369, std::max(j26t, j19t)));
                const int xh = std::min(pm123 + j369, std::min(pm132 + j26t, pm231 + j19t));
                double At = 0;
                for (auto x = xl; x <= xh; ++x)
                {
                    At = -At + unsafe_binomial(x + 1, j123 + 1) * unsafe_binomial(pm123, x - j369) *
                                   unsafe_binomial(pm132, x - j26t) * unsafe_binomial(pm231, x - j19t);
                }
                const int yl = std::max(j456, std::max(j26t, std::max(j258, j48t)));
                const int yh = std::min(pm456 + j26t, std::min(pm465 + j258, pm564 + j48t));
                double Bt = 0;
                for (auto y = yl; y <= yh; ++y)
                {
                    Bt = -Bt + unsafe_binomial(y + 1, j456 + 1) * unsafe_binomial(pm456, y - j26t) *
                                   unsafe_binomial(pm465, y - j258) * unsafe_binomial(pm564, y - j48t);
                }
                const int zl = std::max(j789, std::max(j19t, std::max(j48t, j147)));
                const int zh = std::min(pm789 + j19t, std::min(pm798 + j48t, pm897 + j147));
                double Ct = 0;
                for (auto z = zl; z <= zh; ++z)
                {
                    Ct = -Ct + unsafe_binomial(z + 1, j789 + 1) * unsafe_binomial(pm789, z - j19t) *
                                   unsafe_binomial(pm798, z - j48t) * unsafe_binomial(pm897, z - j147);
                }
                PABC += iphase(xh + yh + zh) * At * Bt * Ct / Pt_de;
            }
            return iphase(dth) * P0 * PABC;
        }

        double norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
        {
            return f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9) *
                   std::sqrt((dj3 + 1.) * (dj6 + 1.) * (dj7 + 1.) * (dj8 + 1.));
        }

        static double lsjj_helper(int l1, int l2, int dj1, int dj2, int J)
        {
            // assume (l1, j1), (l2, j2) is checked in the lsjj function
            const int64_t pj = (dj1 + dj2) / 2;
            const int64_t mj = (dj1 - dj2) / 2;
            if (dj1 > 2 * l1 && dj2 > 2 * l2) // j1 = l1 + 1/2, j2 = l1 + 1/2
            {
                return std::sqrt(double((pj + J + 1) * (pj - J)) / double(2 * dj1 * dj2));
            }
            else if (dj1 > 2 * l1 && dj2 < 2 * l2) // j1 = l1 + 1/2, j2 = l1 - 1/2
            {
                return std::sqrt(double((J - mj + 1) * (J + mj)) / double(2 * dj1 * (dj2 + 2)));
            }
            else if (dj1 < 2 * l1 && dj2 > 2 * l2) // j1 = l1 - 1/2, j2 = l1 + 1/2
            {
                return -std::sqrt(double((J + mj + 1) * (J - mj)) / double(2 * dj2 * (dj1 + 2)));
            }
            else // j1 = l1 - 1/2, j2 = l1 - 1/2
            {
                return std::sqrt(double((pj + J + 2) * (pj - J + 1)) / double(2 * (dj1 + 2) * (dj2 + 2)));
            }
        }

        // S = 0
        static double lsjj_S0(int l1, int l2, int dj1, int dj2, int J) { return lsjj_helper(l1, l2, dj1, dj2, J); }

        // S = 1, J = L
        static double lsjj_S1_0(int l1, int l2, int dj1, int dj2, int J)
        {
            const int64_t pj = (dj1 + dj2) / 2;
            const int64_t mj = (dj1 - dj2) / 2;
            const int64_t pl = l1 + l2;
            const int64_t ml = l1 - l2;
            return (mj * (pj + 1) - ml * (pl + 1)) * lsjj_helper(l1, l2, dj1, dj2, J) / std::sqrt(J * (J + 1));
        }

        // S = 1, J = L + 1
        static double lsjj_S1_p1(int l1, int l2, int dj1, int dj2, int J)
        {
            const int64_t pj = (dj1 + dj2) / 2;
            const int64_t mj = (dj1 - dj2) / 2;
            const int64_t pl = l1 + l2;
            const int64_t ml = l1 - l2;
            const double f0 = J * (2 * J + 1);
            const double fL = (J + mj) * (J - mj) * (J + pj + 1) * (pj - J + 1);
            const double fJ = (J + ml) * (J - ml) * (J + pl + 1) * (pl - J + 1);
            return std::sqrt(fL / f0) * lsjj_helper(l1, l2, dj1, dj2, J - 1) -
                   std::sqrt(fJ / f0) * lsjj_helper(l1, l2, dj1, dj2, J);
        }

        // S = 1, J = L - 1
        static double lsjj_S1_m1(int l1, int l2, int dj1, int dj2, int J)
        {
            const int L = J + 1;
            const int64_t pj = (dj1 + dj2) / 2;
            const int64_t mj = (dj1 - dj2) / 2;
            const int64_t pl = l1 + l2;
            const int64_t ml = l1 - l2;
            const double f0 = L * (2 * L - 1);
            const double fJ = (L + ml) * (L - ml) * (L + pl + 1) * (pl - L + 1);
            const double fL = (L + mj) * (L - mj) * (L + pj + 1) * (pj - L + 1);
            return std::sqrt(fJ / f0) * lsjj_helper(l1, l2, dj1, dj2, J) -
                   std::sqrt(fL / f0) * lsjj_helper(l1, l2, dj1, dj2, L);
        }

        static double lsjj(int l1, int l2, int dj1, int dj2, int L, int S, int J)
        {
            if (dj1 < 1 || dj2 < 1)
                return 0;
            if (std::abs(2 * l1 - dj1) != 1 || std::abs(2 * l2 - dj2) != 1)
                return 0;
            if (!check_couple(2 * l1, 2 * l2, 2 * L))
                return 0;
            if (!check_couple(dj1, dj2, 2 * J))
                return 0;
            if (!check_couple(2 * L, 2 * S, 2 * J))
                return 0;
            if (S == 0)
            {
                return lsjj_S0(l1, l2, dj1, dj2, J);
            }
            else if (S == 1)
            {
                if (L == J)
                    return lsjj_S1_0(l1, l2, dj1, dj2, J);
                else if (J == L + 1)
                    return lsjj_S1_p1(l1, l2, dj1, dj2, J);
                else if (J == L - 1)
                    return lsjj_S1_m1(l1, l2, dj1, dj2, J);
            }
            return 0;
        }

        // Buck et al. Nuc. Phys. A 600 (1996) 387-402
        double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda,
                         double tan_beta = 1.0) const
        {
            // check energy conservation
            const int f1 = 2 * n1 + l1;
            const int f2 = 2 * n2 + l2;
            const int F = 2 * N + L;
            const int f = 2 * n + l;
            if (f1 + f2 != f + F)
                return 0;

            const double sec_beta = std::sqrt(1. + tan_beta * tan_beta);
            const double cos_beta = 1. / sec_beta;
            const double sin_beta = tan_beta / sec_beta;

            const int nl1 = n1 + l1;
            const int nl2 = n2 + l2;
            const int NL = N + L;
            const int nl = n + l;
            const double r1 = unsafe_binomial(2 * nl1 + 1, nl1) / (unsafe_binomial(f1 + 2, n1) * (uint64_t(nl1 + 2) << l1));
            const double r2 = unsafe_binomial(2 * nl2 + 1, nl2) / (unsafe_binomial(f2 + 2, n2) * (uint64_t(nl2 + 2) << l2));
            const double R = unsafe_binomial(2 * NL + 1, NL) / (unsafe_binomial(F + 2, N) * (uint64_t(NL + 2) << L));
            const double r = unsafe_binomial(2 * nl + 1, nl) / (unsafe_binomial(f + 2, n) * (uint64_t(nl + 2) << l));
            const double pre_sum = std::sqrt(r1 * r2 * R * r);
            double sum = 0.;
            for (int fa = 0; fa <= std::min(f1, F); ++fa)
            {
                const int fb = f1 - fa;
                const int fc = F - fa;
                const int fd = f2 - fc;
                if (fd < 0)
                    continue;
                const double t = quick_pow(sin_beta, fa + fd) * quick_pow(cos_beta, fb + fc) *
                                 std::sqrt(unsafe_binomial(f1 + 2, fa + 1) * unsafe_binomial(f2 + 2, fc + 1) *
                                           unsafe_binomial(F + 2, fa + 1) * unsafe_binomial(f + 2, fb + 1));
                for (int la = fa & 0x01; la <= fa; la += 2)
                {
                    const int na = (fa - la) / 2;
                    const int nla = na + la;
                    const double ta =
                        ((uint64_t(2 * la + 1) << la) * unsafe_binomial(fa + 1, na)) / unsafe_binomial(2 * nla + 1, nla);
                    for (int lb = std::abs(l1 - la); lb <= std::min(la + l1, fb); lb += 2)
                    {
                        const int nb = (fb - lb) / 2;
                        const int nlb = nb + lb;
                        const double tb = ((uint64_t(2 * lb + 1) << lb) * unsafe_binomial(fb + 1, nb)) /
                                          unsafe_binomial(2 * nlb + 1, nlb);
                        const int g1 = (la + lb + l1) / 2;
                        const double CGab =
                            unsafe_binomial(g1, l1) * unsafe_binomial(l1, g1 - la) /
                            std::sqrt(unsafe_binomial(2 * g1 + 1, 2 * (g1 - l1)) * unsafe_binomial(2 * l1, 2 * (g1 - la)));
                        for (int lc = std::abs(L - la); lc <= std::min(la + L, fc); lc += 2)
                        {
                            const int nc = (fc - lc) / 2;
                            const int nlc = nc + lc;
                            const double tc = ((uint16_t(2 * lc + 1) << lc) * unsafe_binomial(fc + 1, nc)) /
                                              unsafe_binomial(2 * nlc + 1, nlc);
                            const int G = (la + lc + L) / 2;
                            const double CGac =
                                unsafe_binomial(G, L) * unsafe_binomial(L, G - la) /
                                std::sqrt(unsafe_binomial(2 * G + 1, 2 * (G - L)) * unsafe_binomial(2 * L, 2 * (G - la)));
                            const int ld_min = std::max(std::abs(l2 - lc), std::abs(l - lb));
                            const int ld_max = std::min(fd, std::min(lb + l, lc + l2));
                            for (int ld = ld_min; ld <= ld_max; ld += 2)
                            {
                                const int nd = (fd - ld) / 2;
                                const int nld = nd + ld;
                                const double td = ((uint64_t(2 * ld + 1) << ld) * unsafe_binomial(fd + 1, nd)) /
                                                  unsafe_binomial(2 * nld + 1, nld);
                                const int g2 = (lc + ld + l2) / 2;
                                const double CGcd = unsafe_binomial(g2, l2) * unsafe_binomial(l2, g2 - lc) /
                                                    std::sqrt(unsafe_binomial(2 * g2 + 1, 2 * (g2 - l2)) *
                                                              unsafe_binomial(2 * l2, 2 * (g2 - lc)));
                                const int g = (lb + ld + l) / 2;
                                const double CGbd = unsafe_binomial(g, l) * unsafe_binomial(l, g - lb) /
                                                    std::sqrt(unsafe_binomial(2 * g + 1, 2 * (g - l)) *
                                                              unsafe_binomial(2 * l, 2 * (g - lb)));
                                const int phase = iphase(ld);
                                const double ninej =
                                    f9j(2 * la, 2 * lb, 2 * l1, 2 * lc, 2 * ld, 2 * l2, 2 * L, 2 * l, 2 * lambda);
                                sum += phase * t * ta * tb * tc * td * CGab * CGac * CGbd * CGcd * ninej;
                            }
                        }
                    }
                }
            }
            return pre_sum * sum;
        }

        double dfunc(int dj, int dm1, int dm2, double beta) const
        {
            if (!(check_jm(dj, dm1) && check_jm(dj, dm2)))
                return 0.;
            const int jm1 = (dj - dm1) / 2;
            const int jp1 = (dj + dm1) / 2;
            const int jm2 = (dj - dm2) / 2;
            const int mm = (dm1 + dm2) / 2;
            const double c = std::cos(beta / 2);
            const double s = std::sin(beta / 2);
            const int kmin = std::max(0, -mm);
            const int kmax = std::min(jm1, jm2);
            double sum = 0.;
            for (int k = kmin; k <= kmax; ++k)
            {
                sum = -sum + unsafe_binomial(jm1, k) * unsafe_binomial(jp1, mm + k) * quick_pow(c, mm + 2 * k) *
                                 quick_pow(s, jm1 + jm2 - 2 * k);
            }
            sum = iphase(jm2 + kmax) * sum;
            sum = sum * std::sqrt(unsafe_binomial(dj, jm1) / unsafe_binomial(dj, jm2));
            return sum;
        }

        void reserve(int num, std::string type, int rank)
        {
            if (type == "Jmax")
            {
                switch (rank)
                {
                case 3:
                    fill_binomial_data(3 * num + 1);
                    break;
                case 6:
                    fill_binomial_data(4 * num + 1);
                    break;
                case 9:
                    fill_binomial_data(5 * num + 1);
                    break;
                default:
                {
                    std::cerr << "Error: rank must be 3, 6, or 9" << std::endl;
                    std::exit(-1);
                }
                }
            }
            else if (type == "2bjmax")
            {
                switch (rank)
                {
                case 3:
                    fill_binomial_data(2 * num + 1);
                    break;
                case 6:
                    fill_binomial_data(3 * num + 1);
                    break;
                case 9:
                    fill_binomial_data(4 * num + 1);
                    break;
                default:
                {
                    std::cerr << "Error: rank must be 3, 6, or 9" << std::endl;
                    std::exit(-1);
                }
                }
            }
            else if (type == "nmax")
            {
                fill_binomial_data(num);
            }
            else
            {
                std::cerr << "Error: type must be Jmax, 2bjmax, or nmax" << std::endl;
                std::exit(-1);
            }
        }

    private:
        std::vector<double> _binomial_data;
        int _nmax;
        static std::size_t _binomial_data_size(int n)
        {
            std::size_t x = n / 2 + 1;
            return x * (x + isodd(n));
        }
        static std::size_t _binomial_index(int n, int k)
        {
            std::size_t x = n / 2 + 1;
            return x * (x - iseven(n)) + k;
        }
        void fill_binomial_data(int nmax)
        {
            if (nmax <= _nmax)
                return;
            std::vector<double> old_data = _binomial_data;
            std::size_t reserve_size = _binomial_data_size(nmax);
            if (reserve_size > std::numeric_limits<int>::max())
            {
                std::cerr << "Error: nmax too large" << std::endl;
                std::exit(-1);
            }
            _binomial_data.resize(reserve_size);
            std::copy(std::begin(old_data), std::end(old_data), _binomial_data.begin());
            for (int n = _nmax + 1; n <= nmax; ++n)
            {
                for (int k = 0; k <= n / 2; ++k)
                {
                    _binomial_data[_binomial_index(n, k)] = binomial(n - 1, k) + binomial(n - 1, k - 1);
                }
                ++_nmax;
            }
            _nmax = nmax;
        }
    };

    inline WignerSymbols wigner;

    inline void wigner_init(int num, std::string type, int rank) { wigner.reserve(num, type, rank); }

    inline double fast_binomial(int n, int k) { return wigner.binomial(n, k); }

    // CG coefficient for two spin-1/2
    inline double CGspin(int ds1, int ds2, int S)
    {
        static constexpr double inv_sqrt_2 = 0.70710678118654752;
        if (unsigned(S) > 1)
            return 0;
        if (std::abs(ds1) != 1 || std::abs(ds1) != 1)
            return 0;
        if (S == 0)
        {
            return ds1 == ds2 ? 0.0 : std::copysign(inv_sqrt_2, ds1);
        }
        else // S == 1
        {
            return ds1 == ds2 ? 1.0 : inv_sqrt_2;
        }
    }

    inline double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
    {
        return wigner.CG(dj1, dj2, dj3, dm1, dm2, dm3);
    }

    inline double CG0(int j1, int j2, int j3) { return wigner.CG0(j1, j2, j3); }

    inline double wigner_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
    {
        return wigner.f3j(dj1, dj2, dj3, dm1, dm2, dm3);
    }

    inline double wigner_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
    {
        return wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
    }

    inline double Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
    {
        return wigner.Racah(dj1, dj2, dj3, dj4, dj5, dj6);
    }

    inline double wigner_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
    {
        return wigner.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
    }

    inline double wigner_norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
    {
        return wigner.norm9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
    }

    inline double lsjj(int l1, int l2, int dj1, int dj2, int L, int S, int J)
    {
        return wigner.lsjj(l1, l2, dj1, dj2, L, S, J);
    }

    inline double dfunc(int dj, int dm1, int dm2, double beta) { return wigner.dfunc(dj, dm1, dm2, beta); }

    inline double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda, double tan_beta = 1.0)
    {
        return wigner.Moshinsky(N, L, n, l, n1, l1, n2, l2, lambda, tan_beta);
    }

} // end namespace util

#endif // JSHL_WIGNERSYMBOL_HPP