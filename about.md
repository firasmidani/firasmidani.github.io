---
layout:    page
title:     About
permalink: /about/
---

Influenced by Vida, a minimalist jekyll theme.

- Author: Firas Midani
- Email:  firasmidani@gmail.com
- Github: [https://github.com/firasmidni](https://github.com/firasmidani)


<div class="blurb">
<table>
<col width="78%">
<col width="22%">
<tr> 
   <td align="left" vertical-align="center" style="padding-left:0px">
         <p1 class="blurb">I am a <a href="https://genome.duke.edu/education/CBB">computational biologist</a> studying the ecology of human gut microbiome. 
	   My scientific research leverages high-throughput experiments and machine learning 
		to illustrate how ecosystem structure influences ecosystem functions such as 
		community productivity, community assembly, and colonization resistance. 
	   I am currently a Ph.D. candidate in the <a href="http://el.ladlab.org:8080">Lawrence David lab</a> 
		at <a href="https://en.wikipedia.org/wiki/Duke_University">Duke University</a>. 
       </p1>  
</td>
<td align="right" style="padding-left:0px">
<img src="/assets/midani_head.png" width="120px" class="img-polaroid" opacity="0.5";/>
</td>  
</tr>
</table>
</div>


<div id="home">

  <h2>Blog Posts</h2>

    <ul class="posts">
    {% for post in site.posts %}
	<li>
		<span>{{ post.date | date: "%F" }} &raquo;</span>  
	        <a style="font-size:1em;" href="{{ post.url }}">{{ post.title }}</a>
	</li>
    {% endfor %}
	  
  </ul>
</div>
