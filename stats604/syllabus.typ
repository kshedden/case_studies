#show link: set text(fill: blue)

= Statistical Practice (Statistics 604)

== Instructor

Kerby Shedden (#link("mailto:kshedden@umich.edu"))\
Office: 277 West Hall\
Office hours: Monday 10:30-12, Wednesday 3:30-5

== Overview

This is a PhD level course about the practice of statistics.  It is intended for students in the second year of the Statistics PhD program who have previously completed Statistics 600 and either 511 or 610.

Statistical practice encompasses foundational knowledge and practical skills enabling rigorous data analysis. At a high level, it includes

- Formulating compelling and tractable research questions
- Understanding the capacity of a particular dataset and specific statistical methods to address the research questions
- Implementing statistical methods and conducting the analysis
- Communication of methods and findings

Our course will use a series of case studies to help students build this foundation of knowledge and skills.  Each case study is based on open-access data.  We will discuss each case study for 2-3 weeks, and therefore we will cover around 6 case studies during the term.  These case studies involve complex data and scientific challenges in diverse fields of inquiry.    By having a common set of case studies, we will be able to engage in a shared, focused discussion about a specific set of ideas, both methodological and scientific.

In addition to working with the data, you will be assigned readings throughout the term, which will mostly be papers from the statistics literature, or from the broader scientific literature that involves data analysis.

This is a small, seminar-style class in which participation and attendance are both very important.  The instructor will lead most class discussions, but continuous engagement and frequent participation of all students in these discussions is essential, and will be an important part of the assessment.

== Learning goals

- Students will be able to rapidly internalize (at an appropriately high level), the current state of scientific knowledge in a research domain, enabling them to participate in the formulation of research aims.

- Students will be able to fully assess the structure of a dataset and its capacity for addressing a research aim.  This includes characteristics of the dataset itself as well as all aspects of sampling and measurement that were involved with the production or collection of the data.

- Students will be familiar with a wide variety of methods for data analysis, and deeply understand their appropriate use, strengths, and limitations.  Students will be able to develop an analysis plan that fully justifies how the proposed methods are suited for addressing the scientific research aim.

- Students will be able to assess _a priori_ the potential outcomes and failure modes of a proposed analysis.  This includes an understanding of detectable effect sizes (statistical power), and foreseeable limitations of the research results.

- Students will be adept in communicating the methods and findings of their research to a variety of stakeholders including domain experts, fellow statisticians, persons with expertise in adjacent fields, and the general public.  This will include both written and spoken communication, and will consider both the iterative back-and-forth communication that happens while a project is underway, and how the methods and findings are reported in their final form to various audiences.

== Rationale for the course design

There are many ways that a course like this can be structured.  In the interest of transparency, here is a "Socratic dialogue" about some of the choices made in designing this course.  I also threw in some comments I have heard students make about courses like this in the past and how I typically respond.

- _Why can't we analyze our own data and provide our own research questions?_  If we did this, then everyone would have their own version of the class, and it would be difficult for us to have a discussion about common topics.  Also, as a practicing statistician you are rarely solely empowered to select data and define aims, although it is important to be able to contribute meaningfully to these decisions.

- _I'm not interested in your datasets!_ I admit that the datasets in this course reflect my biases to some extent.  They cover a range of prominent research areas and broadly deal with human-centered topics including health and the environment that affect us all.

- _The methods you discuss are old-fashioned._ We will cover a variety of methods in this course, including some traditional methods like generalized linear models as well as some methods that emerged only in the past 10 years.  When analyzing data in collaboration with domain scientists, we should generally use established methods for the things that they are good at, while thinking creatively about opportunities to reframe the analysis in a way that allows newer methods to contribute.

- _This class has too much/too little theory._ By "theory" here we mean a rigorous assessment of the properties of an analytic method, which almost always involves understanding theoretical results, although in this class we will not be proving any theorems.  This class is more at the applied end of the spectrum, but statistical theory remains highly relevant -- it provides the basis for our understanding of how statistical methods perform, their strengths, and their limitations.

- _Why isn't this more like Kaggle?_ Prediction contests became very popular at the peak of the ML boom (which is passed, as we are now in the AI boom).  While often useful, answering fundamental scientific questions rarely hinges primarily on fitting models that predict well.  Data-driven models that capture important elements of the natural world can be used to answer fundamental scientific questions, but they may not always predict well, and fitting the model is often easier than interpreting and drawing meaningful conclusions from it, which will be the focus of this course.

- _You are too picky about minor details of our writing._ While everyone should write in their own voice, academic writing is also fairly constrained.  When reading student writing I frequently see things that would almost never be found in published academic writing.  Your audience is extremely busy and has many choices of what to read.  If you make your reader work harder than necessary to follow your arguments, they will just give up and find something else to read.  Or, if they are reviewing your work they will evaluate it negatively.

- _I just want to write research papers and this is not helping me do that._ We want to help you launch your career as a scholar of statistics, which may take place in a number of different settings.  The vast majority of our PhD students do an industry internship, and more than half begin their career in industry.  While this course cannot possibly prepare you for all of the different ways that statistics is used in industry, I feel that by working with datasets that address current and deep questions of scientific research, you will gain skills that you can transfer to other settings, including in industry.  

== Coursework and grading

There are three types of work that you will do for this course, and that will determine your grade.

- _Short writing assignments (memos) of approximately two pages_ (50%):  There will be 6-8 of these during the term. For each writing assignment, you will be provided with a prompt and asked to respond to it.  In most cases, this will involve conducting and reporting on analyses of data that we provide you.  But in some cases the prompt may be to comment on a paper that you are assigned to read, or to reflect on some topics that have arisen in class discussions. 

- _Attendance and participation in class discussions_ (15%)

- _In-class writing_ (15%): On October 28th you will write a short memo in class in a blue book (hand-written), responding to a given prompt.  This writing will be done without any reference materials (no notes, books, phones, computers, etc.).

- _Individual oral exams_ (20%): During the last five lab meetings (starting November 4th), each student will have a 15 minute oral exam based on the course materials, tailored to the contents of your individual responses to the writing assignments. The main goal is to assess how well you can explain and defend the contents of your memos, and accurately comment on the major ideas that emerged during class discussions, including in the readings.

== Format

Every course meeting (on Tuesdays and Thursdays) will be a discussion, not a lecture.  Attendance is mandatory and you are expected to be fully engaged in the class discussion.

== Course materials

Course materials will be available through the course
#link("https://github.com/kshedden/case_studies")[github site]. There is no textbook for this course. Announcements and grading will use the course canvas site.

== Case studies and datasets

The case studies will cover a variety of research domains including human health, environment, social, behavioral, physical, and natural science. We will spend approximately three lectures discussing the background and relevant statistical issues associated with each case study. After this context is provided, you will be given a writing prompt and will have approximately 10 days to prepare a report. Extensions will be given only in extreme situations and must be approved in advance of the deadline.

== Computing

You will need to write code to analyze data for this course, but computing is not a major focus of course discussions or evaluation.  Through the course #link("https://github.com/kshedden/case_studies")[github site], we will provide you with a lot of code that can be used to work with the data in each case study.  Depending on your research aims you may need to make only minor modifications and extensions to the code that we provide, but in some cases more extensive coding will be necessary.  You will not be submitting any code to us for grading or evaluation.

== Collaboration, AI usage, and academic honesty

Every memo that you write for this class should reflect your own independent understanding of the ideas that you are presenting and should reflect your own human voice. You alone are responsible for the accuracy and authenticity of your work.  Each student must complete their own analysis and writing, but you may discuss analysis strategies and coding issues with other students. Each student should develop and maintain their own code to conduct the analysis for each case study, although you are welcome to use code provided through the course github site or from other sources as you deem appropriate.

In the submitted writings, plagiarism from the instructor's lectures, course materials, writing of other students, or from any non-acknowledged source is a serious violation of academic integrity and will result at minimum in a score of zero on the assignment, and may also be referred to the LSA or Rackham Honor Council.

You may use large language models to assist in background research, writing, and coding, and you do not need to acknowledge this. However LLMs may produce a generic or superficial analysis of the data that is not reflective of the extended and highly focused discussions that we will be having in class. You need to develop your own understanding of the methods and data that we study in this course, and you need to develop your own distinctly human voice in writing about data analysis. It is strongly suggested that you draft each memo without using an LLM, and then if you wish, you can use the LLM to polish the writing for clarity.

== Expectations for student writing

All statisticians must frequently communicate and document their findings in written form. Writing about your analytic plans and research findings is an excellent way to organize your thoughts, strengthen your arguments, and identify weaknesses in your claims (this is "writing to learn").

The type of writing emphasized in this course reflects the type of writing routinely done by data analysts in many academic and industrial settings. These are short, focused writing assignments and are intended to document a data analysis or provide a concise viewpoint about a scientific or statistical topic. The writing should be rich in specific content. All of the writing should focus on communicating specific scientific insights -- these are not essays written for casual reading by a general audience.

Your writing should focus on a specific question, and should exclusively present analysis and findings that address that question. You generally should provide one short paragraph motivating the research question, and this motivating paragraph may be somewhat less technical than the remainder of the memo. Apart from this brief introduction, your writing should exclusively focus on analysis and findings.  If you are asked to reflect on or critique a paper or idea, you will still want your writing to be focused on a small number of connected main points.

For the writing assignments in this course, you should imagine that you are writing a memo to be read by colleagues or collaborators. Your audience consists of people familiar with the data and scientific or industrial context behind the data, and who are also knowledgeable about statistics at your own level.

Below are some guidelines for writing in this course. These are a few of the most important things to keep in mind. We will expand on this a lot during the semester.

- Write in an academic or professional tone. Do not write casually, but also avoid excessively formal language.

- Organize your content so that your writing has a focused message, with each paragraph contributing a specific and distinct point to this message.

- Each report should focus on one primary question.

- Motivate your analytic approaches in terms of fundamental ideas from
  statistics.

- Use simple, direct language. Avoid convoluted expressions, hidden or
  needlessly subtle meanings, unusual vocabulary, and digressions that do not contribute to your main message.

- Express your arguments in as plain and simple terms as is practical. Do not write in a way that requires the reader to re-read your writing multiple times to understand your point. Favor short sentences in the active voice with a single clause.  Your paragraphs should focus on one topic and generally should be limited to 5-6 sentences.

- Write for an audience that may include non-native English speakers. Avoid colloquialisms and obscure cultural references or metaphors.

- You are expected to be able to write in grammatically-correct English. It is understood that there are many non-native English speakers in the course, but with AI support it is easy to identify and resolve almost all grammatical issues in your writing.