# IBT618-SystemsBiotechnology
-----
### Term project guideline
* All individual projects will be maintained in the [Github `IBT618` repository](https://github.com/igchoi/IBT618-SystemsBiotechnology/2024-spring/)
  * What is `Git` and `Github`?
    * __using Git is not an obligation__
    * [생활코딩 Git](https://opentutorials.org/module/3733/22434)
    * [What is `Github`?](https://www.youtube.com/watch?v=w3jLJU7DT5E)

* Use `markdown` for documentation
  * [Markdown tutorial](https://guides.github.com/features/mastering-markdown/)
  * [Markdown tutorial - Korean](https://github.com/biospin/BigBio/blob/master/reference/%EB%A7%88%ED%81%AC%EB%8B%A4%EC%9A%B4.md)
* You can use [`jupyter notebook`](https://jupyter.org/) too
* Check [library](https://library.korea.ac.kr/) - there are plenty of reference books for python codings
* __Evaluation is on your *activity* rather than your *performance* (e.g. how frequently commit/update your codes/documents/discussions)__ 

### How to install miniconda (to make `conda` environemnt)
* What is [conda]()? 프로그래밍 가상환경
* 설치방법
터미널(terminal)에서 아래 명령어를 입력하고 실행합니다.
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```
* 윈도우즈 사용자는 아래와 같이 리눅스개발환경을 설치하면 편리합니다.
#### 윈도우즈 사용자 - 리눅스개발환경 설치
* 윈도우즈 사용자는 리눅스개발환경을 만드는 것을 권장 - [WSL 설치](https://learn.microsoft.com/ko-kr/windows/wsl/install) 방법에 따라 가상환경을 만듭니다.
```
관리자 모드에서 PowerShell 또는 Windows 명령 프롬프트를 마우스 오른쪽 단추로 클릭하고
"관리자 권한으로 실행"을 선택하여 열고 wsl --install 명령을 입력한 다음 컴퓨터를 다시 시작합니다
```
* WSL 터미널에서 리눅스환경을 이용합니다.


### Links
* [따라하는 데이터과학](https://dataninja.me/ipds-kr/slides-ppt/)
  - `dplyr` package (`tidyverse`)
* [Learn Bioconductor](https://github.com/Bioconductor/LearnBioconductor)

---
### Term projects 
- You can choose following topics or suggest your own term project
- In following weeks, everyone will present and discuss about his/her intro/progress/issues/results, including 
  - i) what is the problem to solve? (_topic description/introduction_), 
  - ii) why is it important? (_biological meaning_), 
  - iii) how it can be solved (_algorithms/procedures_) - _how do you solve the problem?_ 
- Some projects are mostly related to `indexing` and `pattern matching` of sequences, which was handled in the previous class hours.
