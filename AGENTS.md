
# Study Details
This is an fNIRs study investigating the two way interaction of short-form videos impact on attention and retention. The study is organized as a within-subjects study with 4 conditions each repeated 4 times.
The conditions are: short-form education, short-form entertainment, long-form education, and long-form entertainment. We want to evaluate the prefrontal cortex activation for the impact of length and content type on brain activation and again on a participants ability to recall information afterwards.

More information on the study can be found in README.md.

# Scientific Integrity
This study is due to be published in a major cognitive science journal and as such the quality of the preprocessing, statistical analysis, results, etc needs to be of the highest scientific quality. If you implement something please always cite the peer-reviewed source from which you derived the information. Everything must be defensible with citations both immediately in the comment of the section and in a CITATIONS.md file. The CITATIONS.md file needs to be cleanly organized by section (i.e preprocessing, statistical analysis, etc) with the authors name, how the citation was used, a link to the work, and reasoning for the use as it relates to this study.


# Vocabulary
- "Run" indicates the whole fNIRs sequence which is made up of a any number of events or tasks.
- "Task" is a subset of a run the start of which is marked by an event or trigger. Ours in this experiment are "1"..."4".
- "Trials" is synonymous with Task

Do:
- Always ensure the documentation is up to date. We never want out of date .md files that pollute the codebase. 
- Cite your sources
- Write clean and maintainable code with function doc strings and comments that enable quick and easy parsing of the code by a human.
- Prioritize modularity and single responsibility when possible. Modular functions with explicit tests are integral to trusting results.
- Test for scientific quality and integrity vs if a program will crash or not. I prefer the scripts to fail hard and fast than to introduce subtle bugs that impact the downstream results.
    - Scientific quality means that if a function says it's preprocessing the spike artifacts, make sure it's doing exactly that. If a statistical test is intended to test a one-sided t-test use a true positive and a true negative test to ensure the test is performing the operation exactly as intended. Be incredibly thorough with the tests. The deeper they go the better.
- Push back on things if I say or ask for something that goes against the scientific norm. Educate me on why the thing I asked for is not likely the right option.

Do not:
- Make assumptions about the analysis. Always stop and ask.
- Introduce code that can silently change the results (i.e fill NaN with 0 or something like that.)